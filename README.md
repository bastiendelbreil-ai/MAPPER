# M.A.P.P.E.R. 🧬
**Motif Analyzer and Promoter Predictor for Expressions and Regulation**

M.A.P.P.E.R. est une application web locale développée en Python (via Streamlit) destinée à l'analyse de génomes bactériens. Elle permet d'automatiser la recherche d'ORFs (Open Reading Frames) et d'identifier la structure complète des opérons en scannant les séquences promotrices, les sites de liaison aux régulateurs (RBS), et les terminateurs transcrits.

Initialement calibré avec des matrices de poids positionnel (PSSM) par défaut pour *Bacillus*, l'outil est universel et accepte l'importation de matrices personnalisées pour scanner n'importe quel contig bactérien issu d'un séquençage (format FASTA/Multi-FASTA).

## 🚀 Fonctionnalités
* **Recherche Multi-FASTA** : Analyse complète de l'ensemble des contigs d'un fichier de séquençage.
* **Notation de Confiance** : Attribution d'un score dynamique (/100) pour chaque gène et opéron basé sur la qualité des signaux (Boîtes -35/-10, RBS, Biais GC3, Thermo-stabilité du terminateur).
* **Cartographie Visuelle** : Génération de cartes génomiques interactives pour chaque opéron détecté.
* **Matrices Personnalisées dynamiques** : Chargement de séquences de régulateurs, de facteurs Sigma et d'espacements (spacers) spécifiques à votre modèle biologique via un fichier Excel `.xlsx`.
* **Exportation** : Sauvegarde des résultats complets (positions, séquences nucléotidiques et protéiques) en un clic au format Excel.

## 🛠️ Installation et Prérequis

L'application fonctionne entièrement en local. Aucune connexion internet n'est requise après l'installation, garantissant la confidentialité de vos données génomiques. Vous devez disposer de **Python 3.8+** sur votre machine.

1. **Cloner le dépôt :**
  `git clone [https://github.com/votre-nom/MAPPER.git](https://github.com/votre-nom/MAPPER.git)`
  `cd MAPPER`

**Installer les bibliothèques requises :**
Il est recommandé de créer un environnement virtuel (ex: via conda). Installez ensuite les dépendances :
`pip install -r requirements.txt`

## 🖥️ Utilisation
Pour lancer l'application, ouvrez votre terminal dans le dossier contenant le script et exécutez :
`streamlit run sequence_analyzer.py`
*(Remarque : remplacez sequence_analyzer.py par le nom exact de votre fichier Python s'il est différent).*
Une fenêtre s'ouvrira automatiquement dans votre navigateur web par défaut à l'adresse locale localhost:8501.

**Étapes de l'analyse :**
1. Paramétrage : Ajustez la taille minimale des ORFs et les seuils d'affinité pour les promoteurs/régulateurs dans le menu latéral.

2. Matrices (Optionnel) : Si vous travaillez sur une autre bactérie, téléchargez le template Excel en cliquant sur le bouton dédié dans l'interface. Remplissez-le avec vos séquences, puis chargez-le. Si aucun fichier n'est fourni, les matrices par défaut seront utilisées.

3. Scan : Glissez-déposez votre fichier génomique (.fasta, .fa, .fna). L'analyse se lance automatiquement en temps réel (mise en cache pour plus de rapidité).

4. Export : Une fois le calcul terminé, naviguez visuellement contig par contig ou cliquez sur "Télécharger le tableau complet (Excel)" pour extraire toutes vos données.
