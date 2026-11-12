FROM node:22-slim

LABEL maintainer="Noah T. Jones"
LABEL description="rNMP results browser (Authentik-gated). v1 is wild type vs RNase H2."

WORKDIR /app

COPY package.json package.json
RUN npm install --production

COPY server.js server.js
COPY public/ public/

EXPOSE 8050

CMD ["node", "server.js"]
