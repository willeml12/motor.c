# Simulation d'un moteur à reluctance variable par la méthode des éléments finis
Dans le cadre du cours [LEPL1110 - Eléments finis](https://sites.uclouvain.be/archives-portail/cdc2020/cours-2020-lepl1110), les étudiants ont implémenté un programme permettant de simuler un moteur à reluctance variable réprésenté par un maillage (encrypté sous format .txt) via la méthode des éléments finis. Ceci est le version fonctionnelle et partiellement optimisée rendue par Laureline Willem.
## Table de matière
* [Environnement](#environnement)
* [Compilation](#compilation)
* [Exécution](#exécution)
* [Bibliographie](#bibliographie)
* [Auteur](#auteur)
## Répertoire et fichier
*  /build contient le dernier build [Visual Studio](https://visualstudio.microsoft.com/fr/).
*  /deps contient la librairie Basic OpenGL Viewer appartenant à Célestin Marot, sous license zlib, nécessaire à l'éxécution de l'interface graphique.
*  /data contient les données fournies en exemple. C'est dans ce dossier que sera ajouté les données de l'utilisateur le cas échant.
*  /src contient le code source
*  CMakeLists.txt permet de genérer un projet avec CMake.
*  README.md, ce fichier.
*  motor.pdf vise à expliquer différents choix d'implémentation et décrire les tests réalisés afin de valider le code et déterminer l'ordre de précision des calculs.
## Environnement
Le code n'a été testé que sous Windows 10 et Ubuntu, aucun garantie ne peut donc être fourni si un utilisateur souhaite l'utiliser sur un autre OS.
## Compilation
Le code peut être compilé via l'usage des outils [Visual Studio](https://visualstudio.microsoft.com/fr/) et [CMake](https://cmake.org/). Il faut alors suivre les consignes d'utilisation sur CMake dans le dossier d'installation souhaite avec pour cible CMakeList.txt et le dossier build. Il est nécessaire de supprimer tout les fichier contenu dans build au préalable pour utiliser cette méthode. Ce n'est cependant pas l'unique solution, tout import du projet dans [Visual Studio](https://visualstudio.microsoft.com/fr/) ou compilation en ligne de commande fonctionnera également.
## Exécution
Les données, rédigée au format attendu par la fonction motorMesh \*motorMeshRead(const char \*filename) disponible à la ligne 277 du fichier /src/main.c.
Il faudrait alors remplacer le chemins vers les données à la ligne 22 du même fichier, l'utilisateur pourra alors exécuter le programme (par exemple via l'outils prévu à cet effet de [Visual Studio](https://visualstudio.microsoft.com/fr/)) pour obtenir une simulation du moteur donné en entrée.
## Bibliographie
* Canevas du projet fourni par l'équipe enseignante LEPL1110
* [Basic OpenGL Viewer](https://git.immc.ucl.ac.be/marotc/anm/-/tree/master/deps/BOV)
* Syllabus du cours LEPL1110 par Vincent Legat et Jean-François Remacle, édition 2020
* [Valgrind](https://valgrind.org/)
## Auteur
Ce programme de simulation d'un moteur à reluctance variable par la méthode des éléments finis à été réaliser par Laureline Willem dans le cadre de l'unité [LEPL1110 - Eléments finis](https://sites.uclouvain.be/archives-portail/cdc2020/cours-2020-lepl1110) comme projet de fin d'année pendant l'année académique 2020-2021.
