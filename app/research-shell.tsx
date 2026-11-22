import { MobileNavigation } from "./mobile-navigation";
import { docsUrl, repoUrl, sitePath, conceptDoi } from "./site-metadata";
import styles from "./research.module.css";

const links = [
  ["CRISPR counting", sitePath("crispr-guide-counting")],
  ["Sensitivity", sitePath("assignment-sensitivity")],
  ["Library checker", sitePath("tools/library-safety")],
] as const;

export function ResearchHeader() {
  return (
    <>
      <a className="skip-link" href="#main-content">
        Skip to main content
      </a>
      <header className="site-header">
        <div className="header-inner">
          <a className="brand" href={sitePath()} aria-label="DotMatch home">
            <span className="brand-mark" aria-hidden="true" />
            <span>DotMatch</span>
          </a>
          <nav className="desktop-nav" aria-label="Primary navigation">
            {links.map(([label, href]) => (
              <a key={href} href={href}>
                {label}
              </a>
            ))}
            <a href={docsUrl}>Docs</a>
            <a href={repoUrl}>GitHub</a>
          </nav>
          <MobileNavigation links={links} docsUrl={docsUrl} repoUrl={repoUrl} />
        </div>
      </header>
    </>
  );
}

export function ResearchFooter() {
  return (
    <footer className={styles.footer}>
      <p>
        <strong>DotMatch</strong>
        <br />
        Local, open-source short-DNA read assignment and QC (Apache-2.0).
      </p>
      <nav aria-label="Footer navigation">
        <a href={`${repoUrl}/blob/main/LICENSE`}>License</a>
        <a href={conceptDoi}>Cite DotMatch</a>
        <a href={`${repoUrl}/issues/new/choose`}>Support and feedback</a>
        <a href={`${docsUrl}trust-and-scope.html`}>Methods and limits</a>
      </nav>
    </footer>
  );
}
