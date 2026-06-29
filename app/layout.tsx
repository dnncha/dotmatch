import type { Metadata, Viewport } from "next";
import "./globals.css";

const siteUrl = process.env.NEXT_PUBLIC_SITE_URL ?? "https://dnncha.github.io/dotmatch";
const socialImageUrl = `${siteUrl}/dotmatch-og.png`;
const twitterImageUrl = `${siteUrl}/dotmatch-twitter.png`;
const socialImageAlt =
  "DotMatch preview showing assignment reliability outcomes for known-target sequencing assays";

export const metadata: Metadata = {
  metadataBase: new URL(siteUrl),
  applicationName: "DotMatch",
  title: "DotMatch - Assignment Reliability for Known-Target Sequencing Assays",
  description:
    "DotMatch shows which known-target read assignments are unique, ambiguous, unmatched, or invalid.",
  authors: [{ name: "DotMatch maintainers", url: "https://github.com/dnncha/dotmatch" }],
  creator: "DotMatch maintainers",
  publisher: "DotMatch",
  category: "Bioinformatics software",
  alternates: {
    canonical: siteUrl
  },
  keywords: [
    "bioinformatics",
    "computational biology",
    "CRISPR",
    "FASTQ",
    "known-target sequencing",
    "assignment reliability",
    "barcode demultiplexing",
    "barcode panel design",
    "barcode troubleshooting",
    "FASTQ sequence matching",
    "nf-core",
    "MultiQC",
    "Galaxy workflows",
    "Snakemake",
    "core facility sequencing"
  ],
  openGraph: {
    title: "DotMatch - Assignment Reliability",
    description:
      "Know which read assignments you can trust for known-target sequencing assays.",
    type: "website",
    siteName: "DotMatch",
    locale: "en_US",
    url: siteUrl,
    images: [
      {
        url: socialImageUrl,
        secureUrl: socialImageUrl,
        width: 1200,
        height: 630,
        type: "image/png",
        alt: socialImageAlt
      }
    ]
  },
  twitter: {
    card: "summary_large_image",
    title: "DotMatch - Assignment Reliability",
    description:
      "Know which read assignments you can trust for known-target sequencing assays.",
    images: [
      {
        url: twitterImageUrl,
        width: 1200,
        height: 630,
        alt: socialImageAlt
      }
    ]
  }
};

export const viewport: Viewport = {
  width: "device-width",
  initialScale: 1
};

export default function RootLayout({
  children
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  );
}
