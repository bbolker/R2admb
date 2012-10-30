(TeX-add-style-hook "R2admb"
 (lambda ()
    (LaTeX-add-bibliographies)
    (LaTeX-add-labels
     "fig:rfsp1")
    (TeX-add-symbols
     '("code" 1)
     '("fixme" 1)
     "R"
     "Splus"
     "windows")
    (TeX-run-style-hooks
     "natbib"
     "hyperref"
     "verbatim"
     "listings"
     "fancyvrb"
     "alltt"
     "url"
     "inputenc"
     "utf8"
     "babel"
     "american"
     "latex2e"
     "art11"
     "article"
     "11pt")))

