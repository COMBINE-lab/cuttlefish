// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import remarkGfm from 'remark-gfm';

// GitHub Pages: site + base. The repo deploys to
// https://combine-lab.github.io/cuttlefish/.
export default defineConfig({
  site: 'https://combine-lab.github.io',
  base: '/cuttlefish',
  // GitHub-Flavored Markdown (tables, strikethrough, autolinks) is applied to
  // `.md` by default but NOT to `.mdx`; add remark-gfm explicitly so tables in
  // the `.mdx` pages render.
  markdown: {
    remarkPlugins: [remarkGfm],
  },
  integrations: [
    starlight({
      title: 'Cuttlefish',
      description:
        'Cuttlefish — fast, parallel, low-memory construction of compacted de Bruijn graphs from reference sequences or sequencing reads.',
      logo: {
        light: './src/assets/cuttlefish-logo.svg',
        dark: './src/assets/cuttlefish-logo-dark.svg',
        alt: 'Cuttlefish',
      },
      social: [
        {
          icon: 'github',
          label: 'GitHub',
          href: 'https://github.com/COMBINE-lab/cuttlefish',
        },
      ],
      editLink: {
        baseUrl: 'https://github.com/COMBINE-lab/cuttlefish/edit/rust-rewrite/website/',
      },
      // Two doorways, mirroring the split on the landing page: the Rust
      // Cuttlefish 3, which is where Cuttlefish 3 is heading, and the C++ line
      // (Cuttlefish 1 and 2, which is what most published work used, plus the
      // C++ Cuttlefish 3 developed in parallel with the Rust one).
      sidebar: [
        {
          label: 'Cuttlefish 3 (Rust)',
          items: [
            { label: 'Introduction', slug: 'rust/introduction' },
            { label: 'Installation', slug: 'rust/installation' },
            { label: 'Quick start', slug: 'rust/quick-start' },
            { label: 'Command line', slug: 'rust/command-line' },
            { label: 'Output formats', slug: 'rust/output' },
            { label: 'Comparing graphs', slug: 'rust/comparing-graphs' },
            { label: 'Resource control', slug: 'rust/resources' },
            { label: 'Architecture', slug: 'rust/architecture' },
          ],
        },
        {
          label: 'Cuttlefish 1 & 2 (C++)',
          items: [
            { label: 'Overview', slug: 'cpp/overview' },
            { label: 'Installation', slug: 'cpp/installation' },
            { label: 'Usage', slug: 'cpp/usage' },
            { label: 'Output formats', slug: 'cpp/output-formats' },
            { label: 'Example usage', slug: 'cpp/examples' },
            { label: 'Larger k-mer sizes', slug: 'cpp/large-k' },
            { label: 'Differences between 1 & 2', slug: 'cpp/differences' },
            { label: 'The C++ Cuttlefish 3', slug: 'cpp/cuttlefish3' },
          ],
        },
        {
          label: 'About',
          items: [
            { label: 'Which version?', slug: 'about/which-version' },
            { label: 'Citations', slug: 'about/citations' },
            { label: 'Licenses', slug: 'about/licenses' },
          ],
        },
      ],
    }),
  ],
});
