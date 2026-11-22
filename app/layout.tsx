import type { Metadata, Viewport } from "next";
import { JetBrains_Mono, Plus_Jakarta_Sans } from "next/font/google";
import { pageMetadata, siteUrl, repoUrl } from "./site-metadata";
import "./globals.css";

const sans = Plus_Jakarta_Sans({
  subsets: ["latin"],
  weight: ["400", "500", "600", "700", "800"],
  variable: "--font-sans",
  display: "swap",
});

const mono = JetBrains_Mono({
  subsets: ["latin"],
  weight: ["400", "500", "600"],
  variable: "--font-mono",
  display: "swap",
});

export const metadata: Metadata = {
  ...pageMetadata(
    "DotMatch: CRISPR guide counting and barcode QC",
    "Assign FASTQ reads to known CRISPR guides and inline barcodes locally. Export MAGeCK-compatible counts while keeping ambiguous, unmatched, and invalid reads visible.",
  ),
  metadataBase: new URL(siteUrl),
  applicationName: "DotMatch",
  authors: [{ name: "DotMatch maintainers", url: repoUrl }],
  creator: "DotMatch maintainers",
  publisher: "DotMatch",
  category: "Bioinformatics software",
  icons: { icon: [{ url: `${siteUrl}/favicon.svg`, type: "image/svg+xml" }] },
};

export const viewport: Viewport = { width: "device-width", initialScale: 1 };

export default function RootLayout({
  children,
}: Readonly<{ children: React.ReactNode }>) {
  return (
    <html lang="en" className={`${sans.variable} ${mono.variable}`}>
      <head>
        <link rel="describedby" href={`${siteUrl}/llms.txt`} />
      </head>
      <body>{children}</body>
    </html>
  );
}
