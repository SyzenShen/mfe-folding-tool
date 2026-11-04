FROM python:3.11-slim

# Non-root user with uid 1000 (also what Hugging Face Docker Spaces expect).
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH \
    PYTHONUNBUFFERED=1 \
    MFE_MAX_SEQS=20 \
    MFE_MAX_LEN=200
WORKDIR $HOME/app

COPY --chown=user requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt gunicorn

COPY --chown=user app.py .

# Hosting platforms set $PORT (Render, Cloud Run, ...); 7860 is the Hugging Face default.
# 2 worker processes x 2 threads; ViennaRNA global state is guarded by a lock in app.py.
CMD ["sh", "-c", "exec gunicorn -b 0.0.0.0:${PORT:-7860} --workers 2 --threads 2 --timeout 120 app:app"]
