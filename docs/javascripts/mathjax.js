window.MathJax = {
  loader: {
    load: ["[tex]/boldsymbol"],
  },
  tex: {
    packages: { "[+]": ["boldsymbol"] },
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true,
    tags: "ams",
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex",
  },
  startup: {
    typeset: false,
  },
};

document$.subscribe(() => {
  MathJax.startup.promise.then(() => {
    MathJax.startup.output.clearCache();
    MathJax.typesetClear();
    MathJax.texReset();
    return MathJax.typesetPromise();
  });
});
