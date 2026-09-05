import type { NextConfig } from "next";
const basePath = process.env.NEXT_PUBLIC_BASE_PATH || undefined;
const nextConfig: NextConfig = {
  output: process.env.NEXT_OUTPUT === "export" ? "export" : undefined,
  // GitHub Pages needs directory/index.html routes for canonical trailing-slash URLs.
  trailingSlash: true,
  basePath, assetPrefix: basePath, devIndicators: false,
  images: { unoptimized: true }
};
export default nextConfig;
