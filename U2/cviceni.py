# proženeme algoritmus, který nám získá korelace pro jednotlivé c, r -> bude to vytvářet shluky (jedna obce s kostelem bude x-krát, ale někde bude mít největší korelace)
# najdeme 5 vzorů s kostelem -> průměrný vzor -> nikde nebude kor. koef. 1, ale bude podobnější průměrné obci
# vypočítat korelační koeficient (O, 1), který se ukáže na testovacím obrázku
# vytvořit ty obdelníčky na původní mapě, ohraničující obce s kostelem

# načíst soubor
# normalizace RGB složek na (0, 1) [rgb/256]

# proměnné pro widht, height vzoru
# normxcorr2! velice důležitý (vstup: vzor, celý obrázek), lze pouze dát jednu barevnou složku

# nastavení hodnoty korelace (cca 0,6)
# for cyklus s porovnáváním vypočtených hodnot korelace s zadaným práhem

# následné vykreslení obdelníčku

# následně tedy vybrat polohu, kde je v jednotlivých shlucích největší kor. koef.

# nějakej program irfan pro odečtení pixelu? (já použiju kritu)

# vhodné použít HSL (L)

# skupina č.2