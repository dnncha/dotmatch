import type { Metadata, Viewport } from "next";
import { pageMetadata, siteUrl, repoUrl } from "./site-metadata";
import "./globals.css";
export const metadata: Metadata = {
  ...pageMetadata("DotMatch — CRISPR guide counting & barcode QC", "Count CRISPR guides from FASTQ, check barcode collisions, and keep ambiguous reads visible. Local, open-source assignment to known short DNA targets."),
  metadataBase: new URL(siteUrl), applicationName: "DotMatch",
  authors: [{ name: "DotMatch maintainers", url: repoUrl }], creator: "DotMatch maintainers", publisher: "DotMatch", category: "Bioinformatics software",
  icons: { icon: [{ url: `${siteUrl}/favicon.svg`, type: "image/svg+xml" }] }
};
export const viewport: Viewport = { width: "device-width", initialScale: 1 };
export default function RootLayout({ children }: Readonly<{ children: React.ReactNode }>) {
  return <html lang="en"><head><link rel="describedby" href={`${siteUrl}/llms.txt`} /></head><body>{children}</body></html>;
}
