import type { Metadata } from "next";
export const siteUrl = (process.env.NEXT_PUBLIC_SITE_URL ?? "https://dnncha.github.io/dotmatch").replace(/\/+$/, "");
export const basePath = (process.env.NEXT_PUBLIC_BASE_PATH ?? "").replace(/\/+$/, "");
export const docsUrl = "https://dotmatch.readthedocs.io/en/latest/";
export const repoUrl = "https://github.com/dnncha/dotmatch";
export const publishedVersion = "0.4.1";
export const conceptDoi = "https://doi.org/10.5281/zenodo.20541628";
export function sitePath(path = ""): string {
  const clean = path.replace(/^\/+|\/+$/g, "");
  return `${basePath}/${clean ? `${clean}/` : ""}`;
}
export function canonicalUrl(path = ""): string {
  const clean = path.replace(/^\/+|\/+$/g, "");
  return `${siteUrl}/${clean ? `${clean}/` : ""}`;
}
export function pageMetadata(title: string, description: string, path = ""): Metadata {
  return {
    title, description, alternates: { canonical: canonicalUrl(path) },
    openGraph: { title, description, type: "website", siteName: "DotMatch", url: canonicalUrl(path),
      images: [{ url: `${siteUrl}/dotmatch-og.png`, width: 1200, height: 630, alt: "DotMatch: known-target read assignment with visible ambiguity and QC" }] },
    twitter: { card: "summary_large_image", title, description, images: [`${siteUrl}/dotmatch-twitter.png`] }
  };
}
