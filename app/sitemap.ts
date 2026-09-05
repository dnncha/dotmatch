import type { MetadataRoute } from "next";
import { canonicalUrl } from "./site-metadata";
export const dynamic = "force-static";
export default function sitemap(): MetadataRoute.Sitemap {
  // Do not manufacture lastModified from build time. Add actual content dates when tracked.
  return ["", "crispr-guide-counting", "tools/library-safety", "assignment-sensitivity"].map(path => ({ url: canonicalUrl(path) }));
}
