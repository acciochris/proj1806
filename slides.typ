#import "@preview/touying:0.7.3": *
#import themes.university: *

#import "@preview/numbly:0.1.0": numbly

#show: university-theme.with(
  aspect-ratio: "16-9",
  config-info(
    title: "Finding Eigenvalues with\nthe QR Algorithm",
    subtitle: "ES.1806 Final Project",
    author: "Chris Liu",
    date: datetime.today(),
    institution: "MIT",
  ),
)

#set heading(numbering: numbly("{1}.", default: "1.1"))

#title-slide()
