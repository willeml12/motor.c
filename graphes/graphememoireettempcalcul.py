import matplotlib.pyplot as plt

NbrNoeudFull = [400, 838, 1667, 4424, 14608]
TempBandY = [12.3, 39.9, 125.5, 802.2, 8505.1]
MemBandY = [69, 70, 72, 81, 143]
TempBandX = [12.5, 44.4, 139.7, 970.6, 9158.6]
MemBandX = [69, 69, 72, 81, 147]
TempBasic = [75.4, 657.2,5961.3]
MemBasic = [72, 75, 92]
NbrNoeudRed = [400, 838, 1667]

plt.title('Temps de calcul moyen en fonction du nombre de noeuds')
plt.xlabel('Nombres de noeuds[#]')
plt.ylabel('Temps[ms]')
plt.plot(NbrNoeudFull, TempBandY, color='teal')
plt.plot(NbrNoeudFull, TempBandX, color='darkorchid')
plt.plot(NbrNoeudRed, TempBasic, color='crimson')
plt.legend(['Solveur Bande (FEM_YNUM)', 'Solveur Bande (FEM_XNUM)', 'Solveur Basique'])
plt.show()

plt.title('Mémoire maximale utilisée en fonction du nombre de noeuds')
plt.xlabel('Nombres de noeuds[#]')
plt.ylabel('Mémoire maximale utilisée[Mo]')
plt.plot(NbrNoeudFull, MemBandY, color='teal')
plt.plot(NbrNoeudFull, MemBandX, color='darkorchid')
plt.plot(NbrNoeudRed, MemBasic, color='crimson')
plt.legend(['Solveur Bande', 'Solveur Bande (FEM_XNUM)', 'Solveur Basique'])
plt.show()