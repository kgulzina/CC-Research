\NeedsTeXFormat{LaTeX2e} 
\ProvidesClass{ejasa}[2013/10/02 LaTeX style for EJASA -- Electronic Journal of Applied Statistical Analysis, v. 0.2]

\LoadClass[11pt,headsepline=false,twoside=semi]{scrartcl}

\RequirePackage{latexsym}
\RequirePackage{bm}
\RequirePackage{multirow}
\RequirePackage{url,amsbsy,amsfonts,amsmath,amssymb,enumerate,epsfig,graphicx,rotating,subfigure}
\RequirePackage[authoryear]{natbib}
\RequirePackage{fancyhdr}
\renewcommand{\arraystretch}{1.2}
\RequirePackage{float}
\RequirePackage[ ]{authblk}
\RequirePackage{booktabs}
\RequirePackage{caption}

\renewcommand\Affilfont{\itshape\small}
%\renewcommand\Affilfont{\itshape}
\addtokomafont{disposition}{\rmfamily}

\newcommand{\covermat}[1]{{\usekomafont{title} #1}}
\newcommand{\jname}{Electronic Journal of Applied Statistical Analysis}

%%%Bibliography: set the 
\setlength\bibsep{0.5\itemsep}

%%% Herader / Footer styles
\setlength{\headheight}{30.0pt}


\fancypagestyle{firstpage}{%
%	\pagenumbering{gobble}
	\fancyhead{}
	\fancyhead[L]{{\itshape{\jname}}\\%
	Vol. \ejasavolume, Issue \ejasanumber, \ejasayear, \ejasafirstpage-\ejasalastpage\\
         DOI: \ejasadoi}
        \fancyfoot[L]{\copyright Universit\`{a} del Salento\\%
        ISSN:  2070-5948\\%
        http://siba-ese.unisalento.it/index.php/ejasa/index
    }

}


\pagestyle{fancy}               %% ...
\fancyhead{}                    %% initialisieren der Kopfzeile
\fancyhead[CE]{\itshape{{\authorlastnames}}}
\fancyhead[LE]{\thepage}
\fancyhead[CO]{\itshape{{\jname}}}
\fancyhead[RO]{\thepage}
\renewcommand{\headrulewidth}{0pt}
%%============ Autorenname

%\fancyhead[LO]{\sffamily{\large{\textcolor{black}{\shortarticletitle}}}}
%%============ Dokumenttitel

\fancyfoot{}                  %% initialisieren der Fu{\ss}zeile
%\fancyfoot[RE]{\sffamily{\bfseries{\textcolor{black}
%{JLCL}}}}  %%<========== links/innen LDV...
%
%%\fancyfoot[LO]{\sffamily{\bfseries{\textcolor{black}
%%{JLCL \ldvforumyear ~-- Band \ldvforumvolume$\,$(\ldvforumnumber)}}}}        %%<========== rechts/innen
%
%\fancyfoot[LE]{\sffamily{\bfseries{\newstylenums{\LARGE{\textcolor{black}
%{\thepage}}}}}}     %%<========== Seitenzahl
%
%\fancyfoot[RO]{\sffamily{\bfseries{\newstylenums{\LARGE{\textcolor{black}
%{\thepage}}}}}}     %%<========== Seitenzahl
%

\endinput
