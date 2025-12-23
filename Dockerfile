# 1. Use the official pixi base image
# 'latest' is usually based on Ubuntu (currently noble or jammy)
FROM ghcr.io/prefix-dev/pixi:latest

# 2. Set the working directory
WORKDIR /app

# 3. Copy dependencies first (better caching)
COPY pyproject.toml pixi.lock ./

# 4. Install dependencies
# --frozen ensures we stick exactly to the lockfile
RUN pixi install --frozen

# 5. Create an entrypoint script to activate the environment
# The 'pixi shell-hook' command generates the shell code needed to activate the environment.
# We wrap this in a script so every command run in this container uses the pixi python.
RUN pixi shell-hook --shell bash > /shell-hook.sh \
    && echo 'exec "$@"' >> /shell-hook.sh \
    && chmod +x /shell-hook.sh

# 6. Uncomment if you need to copy the rest of the project code, or modify as needed.
# COPY . .

# 7. Set the entrypoint
# This script ensures the environment is active for any command (python, bash, etc.)
ENTRYPOINT ["/bin/bash", "/shell-hook.sh"]

# Default command
CMD ["ipython"]