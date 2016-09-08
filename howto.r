# Neues Erstellen des Paketes [aufwendig, und eigentlich nicht notwendig]
# Zur Sicherheit immer eine alte Kopie aufbewahren

#Für Pakete mit C++ Quellcode einfachste Möglichkeit
library(Rcpp)

# Im Workspace dürfen nur die entsprechenden SMOTE und knn Funktionen vorhanden sein
Rcpp.package.skeleton("SMOTE", list = c(ls()),example_code = F)

# Kopieren des C++ Quellcodes in "src" Ordner
# Kopieren der versteckten ( .abc) Funktionen in die entsprechenden *.r Dateien
# Ausführen des folgenden Befehls im entsprechenden Verzeichnis (hier: Verzeichnis dieser Datei)
compileAttributes("SMOTE")
# Ab hier "normal" weiter

# Bearbeitung der bereits vorhandenen Package Dateien
# *.r und *.rd files entsprechend anpassen, description Datei versionsnummer inkrementieren.
# Weitere Anpassungen nach Bedarf (z.B. bei Umbenennungen von Funktionsnamen muss auch die Dokumentation und der jeweilige Aufruf angepasst werden)

# Auf Kommandozeile im entsprechenden Verzeichnis
# R CMD check SMOTE (1 WARNING, 2 NOTES) [passt]
# R CMD build SMOTE -->erstellt SMOTE_X.Y.tar.gz, wobei X.Y die aktuelle Versionsnummer nach der description Datei ist

# installieren
# Zur Sicherheit: Paket löschen, R neustarten, dann installieren 
install.packages("SMOTE_0.2.tar.gz", repos = NULL, type = "source")
library(SMOTE)


?SMOTE
?lnsmote.separated
?knnsearch

#generate data with 2 features
#majority class: 250 samples
data = replicate(2, runif(250))
#first minority class: 50 samples
minority.class1 = replicate(2, runif(50))
#second minority class: 50 samples
minority.class2 = replicate(2, runif(50))

data = rbind(data, minority.class2)


#compute synthetic data for minority.class1
synthetic = smote.separated(minority = minority.class1, k = 5, o = 4, method = 'e')

synthetic = safelevelsmote.separated(data = data, minority = minority.class1, k = 5, method = 'e')

synthetic = lnsmote.separated(data = data, minority = minority.class1, k = 5, o = 4, method = 'e')

?wrapper

?devtools::check
install.packages("devtools")


?check
?devtools::build("SMOTE")
devtools::document("SMOTE")
devtools::check("SMOTE", cleanup = F)
?devtools::check

?document

inst
