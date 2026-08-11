FROM node:22-slim

LABEL maintainer="Noah T. Jones"
LABEL description="HU-rNMP results browser (Authentik-gated)"

WORKDIR /app

COPY package.json package.json
RUN npm install --production

COPY server.js server.js
COPY public/ public/

EXPOSE 8050

CMD ["node", "server.js"]
