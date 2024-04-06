#include "tools.h"
#include "heap.h" // il faut aussi votre code pour heap.c


// Une fonction de type "heuristic" est une fonction h() qui renvoie
// une distance (double) entre une position de départ et une position
// de fin de la grille. La fonction pourrait aussi dépendre de la
// grille (comme par exemple du nombre de murs rencontrés par le
// segment départ-fin), mais on n'utilisera pas forcément ce dernier
// paramètre. Vous pouvez définir votre propre heuristique.
typedef double (*heuristic)(position,position,grid*);


// Heuristique "nulle" pour Dijkstra.
double h0(position s, position t, grid *G){
  return 0.0;
}


// Heuristique "vol d'oiseau" pour A*.
double hvo(position s, position t, grid *G){
  return fmax(abs(t.x-s.x),abs(t.y-s.y));
}


// Heuristique "alpha x vol d'oiseau" pour A*.
static double alpha=0; // 0 = h0, 1 = hvo, 2 = approximation ...
double halpha(position s, position t, grid *G) {
  return alpha*hvo(s,t,G);
}


// Structure "noeud" pour le tas min Q.
typedef struct node {
  position pos;        // position (.x,.y) d'un noeud u
  double cost;         // coût[u]
  double score;        // score[u] = coût[u] + h(u,end)
  struct node* parent; // parent[u] = pointeur vers le père, NULL pour start
} *node;


// Les arêtes, connectant les 8 cases voisines de la grille, sont
// valuées seulement par certaines valeurs correspondant aux
// différentes textures possibles des cases (il y en a 7). Le poids de
// l'arête u->v, noté ω(u,v) dans le cours, entre deux cases u et v
// voisines est déterminé par la texture de la case finale v. Plus
// précisément, si la case v de la grille est de texture t, le poids
// de u->v vaudra ω(u,v) = weight[t] dont les valeurs numériques
// exactes sont définies ci-après. Notez bien que t est un indice
// ('int'), alors que weight[t] est un coût ('double').
//
// La liste des textutres possibles d'une case est donnée dans
// "tools.h": TX_FREE, TX_WALL, TX_WATER, ... Remarquez que
// weight[TX_WALL]<0 ce qui n'est pas a priori une valuation correcte.
// En effet A* ne marche qu'avec des poids positifs ! Mais ce n'est
// pas un problème, puisqu'en position (x,y), si G.texture[x][y] =
// TX_WALL, alors c'est que le sommet à cette position n'existe pas
// dans le graphe ! Et donc aucune arête ne peut donc être incidente à
// (x,y).

double weight[]={
    1.0,  // TX_FREE
  -99.9,  // TX_WALL
    3.0,  // TX_SAND
    9.0,  // TX_WATER
    2.3,  // TX_MUD
    1.5,  // TX_GRASS
    0.1,  // TX_TUNNEL
};


// Que doit faire la fonction A_star(G,h) ?
//------------------------------------------
//
// Votre fonction A_star(G,h) doit construire un chemin dans la grille
// G, entre la position G.start (s) et G.end (t), selon l'heuristique
// h(). Le chemin doit être calculé selon l'algorithme A* vu en cours
// (utilisez les notes de cours !). L'heuristique h() est une fonction
// à choisir lors de l'appel dans le main() parmi h0(), hvo(),
// halpha(). Vous pouvez aussi définir votre propre heuristique et la
// tester.
//
//
// Que doit renvoyer la fonction A_star(G,h) ?
//---------------------------------------------
//
// Si le chemin n'a pas été trouvé (par exemple si la destination est
// enfermée entre 4 murs ou si G.end est sur un mur), il faut renvoyer
// une valeur < 0.
//
// Sinon, il faut renvoyer le coût du chemin trouvé et remplir le
// champs .mark de G pour que le chemin trouvé puisse être visualisé
// par drawGrid(G) (plutard dans le main). Il faut, par convention, la
// case (x,y) appartient au chemin trouvé ssi G.mark[x][y] = MK_PATH.
// Pour que la visualisation soit complète, faite attention à ce que
// G.end et G.start soient bien marqués comme appartenant au chemin.
//
// Utilisez les touches a,z,+,-,p,c pour gérer la vitesse d'affichage
// et de progression de l'algorithme par exemple. Repportez-vous à
// "tools.h" pour avoir la liste des differentes touches et leurs
// actions, ainsi que les différents marquages possibles G.mark[x][y]
// pour une case (x,y).
//
//
// Comment gérer les ensembles P et Q ?
//--------------------------------------
//
// Pour gérer l'ensemble P, servez-vous du champs G.mark[x][y]. On
// aura G.mark[x][y] = MK_USED ssi le noeud en position (x,y) est dans
// P. Par défaut, ce champs est automatiquement initialisé partout à
// MK_NULL par toute fonction d'initialisation de la grille
// (initGridXXX()).
//
// Pour gérer l'ensemble Q, vous devez utiliser un tas min de noeuds
// (type node) avec une fonction de comparaison (à créer) qui dépend
// du champs .score des noeuds. Pour la fonction de comparaison
// inspirez vous de test_heap.c et faites attention au fait que
// l'expression '2.1 - 2.2' une fois castée en 'int' n'est pas
// négative, mais nulle !
//
// Vous devez utilisez la gestion paresseuse du tas (cf. le paragraphe
// du cours à ce sujet, dans l'implémentation de Dijkstra). Pensez
// qu'avec cette gestion paresseuse, la taille de Q est au plus la
// somme des degrés des sommets dans la grille. On peut être plus fin
// en remarquant que tout sommet, sauf s, aura toujours un voisin au
// moins déjà dans P (son parent). Pour visualiser un noeud de
// coordonnées (x,y) qui passe dans le tas Q vous pourrez mettre
// G.mark[x][y] = MK_FRONT au moment où vous l'ajoutez.
//
// Attention ! Même si cela est tentant, il ne faut pas utiliser la
// marque MK_FRONT pour savoir si un sommet (x,y) se trouve déjà dans
// Q. Cela ne sert pas à grand chose, car s'il est dans Q, vous ne
// pourrez pas déplacer le noeud correspondant pour le mettre au bon
// endroit dans Q en fonction de la mise à jour de son score.


int fcmp_noeud(const void *x, const void *y) {
  const double a = (*(node*)x)->score;
  const double b = (*(node*)y)->score;
  return (a<b)? -1 : (a>b); 
}


node createNode(position pos, double cout, node parent , position t, grid G, heuristic h){
  node noeud = malloc(sizeof(*noeud));
  //node noeud;
  noeud -> pos = pos;
  noeud -> cost = cout;        
  noeud -> score = cout + h(pos,t,&G); 
  noeud -> parent = parent;
  return noeud;
}

bool nodeInHeap(heap h, node V){
  for(int i = 1; i<= h->n; i++) if((V->pos.x == (*(node *) h->array[i])->pos.x ) && (V->pos.y == (*(node *) h->array[i])->pos.y ) ) return true;
  return false;
}

int getPos(heap h, node V){
  if(!heap_empty(h)){
    for(int i = 1; i<= h->n; i++) if((V->pos.x == (*(node *) h->array[i])->pos.x ) && (V->pos.y == (*(node *) h->array[i])->pos.y ) ) return i;
  }
  return -1;
}

double A_star(grid G, heuristic h){

  // Pensez à dessiner la grille avec drawGrid(G) à chaque fois que
  // possible, pour visualiser le comportement de votre algorithme.
  // Par exemple, dès vous ajoutez un sommet à P mais aussi lorsque
  // vous reconstruisez le chemin à la fin de la fonction (marquer les
  // sommets à MK_PATH). Lorsqu'un sommet passe dans Q vous pourrez le
  // marquer à MK_FRONT (dans son champs .mark) pour le distinguer à
  // l'affichage des sommets de P (couleur différente).

  // Après avoir extrait un noeud de Q, il ne faut pas le détruire
  // (free), sous peine de ne plus pouvoir reconstruire le chemin
  // trouvé ! Une fonction createNode() peut simplifier votre code.

  // Les bords de la grille sont toujours constitués de murs (texture
  // TX_WALL) ce qui évite d'avoir à tester la validité des indices
  // des positions (sentinelle). Dit autrement, un chemin ne peut pas
  // s'échapper de la grille.

  // Lorsque que vous ajoutez un élément au tas, pensez à tester la
  // valeur de retour heap_add() afin de détecter d'éventuellement
  // dépassement de capacité du tas, situation qui n'est pas censé
  // arriver. Si cela se produit, l'erreur vient soit du tas qui n'est
  // pas assez grand à son initialisation ou que l'algorithme visite
  // anormalement trop de sommets. De même, lorsque vous extrayer un
  // élément, vérifiez qu'il n'est pas NULL, ce qui n'est pas censé
  // arriver. C'est autant de "segmentation fault" que vous pouvez
  // éviter.


  position s = G.start;
  position t = G.end; 
  
  
  /*Premiere Partie*/
  /*P := ∅*/
  // Par defaut le objet G.Mark[x][y] = MK_NULL par tous les points dans la grid.

  /* Q = {S} */
  //printf("works? %d", G.mark[G.X-1][G.Y]); -> ca montre que il y a seulement G.X-1 * G.Y possibles
  int n = (G.X-1) * G.Y; //
  printf("Dimenstion:(%d,%d) | n = %d\n", G.X, G.Y, n);
  heap Q = heap_create(n,fcmp_noeud); 
  node S = createNode(s, 0, NULL, t, G, h); /*coût[s] := 0 parent[s] := ⊥, score[s] := coût[s] + h(s, t)*/
  printf("S:(%d ,%d) \n", S->pos.x,S->pos.y);
  heap_add(Q,&S);
  G.mark[S->pos.x][S->pos.y] = MK_FRONT;
  //Deuxieme Partie
  int compteur = 1;
  while(!heap_empty(Q)){
    //A)

    void *u =  heap_pop(Q);
    //B)
    node *U = ((node*)u);
    position uPos = (*U)->pos;
    printf("U = (%d,%d)\n", uPos.x,uPos.y);
    printf("t = (%d,%d)\n", t.x,t.y);
    //printf("X avec pointeur: %d\n",(*(node*)u)->pos.x); //Pourqoui ((node*)u)[0] <<revien>> le objet node ? il a relation avec le poiteur ? avec void ? ou le deux ? 
    //En verite les objet d'heap ou tas minimum son de type void, donc la facion de l'appeller c'est (<datatype>*)u pour avoir la direction en memoire et .
    //*(<datatype>*)u pour le valeur.

    if((uPos.x == t.x) && (uPos.y == t.y)){
        return (*U)->score;
    }
    //C)
    G.mark[uPos.x][uPos.y] = MK_USED;
    
    //D)
    //Relire l'algorithme, reviser qui sont les <<etages>> suivantes. 

    //Pour tout voisin v pas en P de u

    printf("Boucle numero: %d\n", compteur);
   
    for(int x = uPos.x-1; x <= uPos.x+1; x++){
      for(int y = uPos.y-1; y <= uPos.y+1; y++){
          if(x >= G.X || y > G.Y) break;
          if(x < 0) x = abs(x);
          if(y < 0) y = abs(y);

          printf("Point: (%d,%d)\n",x,y);
          if(x == uPos.x && y == uPos.y) continue;
          printf("Point in P: %d \n",G.mark[x][y] ==  MK_USED);

          if ( G.mark[x][y] != MK_USED ){
            position v = {x,y};
            double cout = sqrt(pow(v.x-uPos.x,2) + pow(v.y-uPos.y,2));
            node V = createNode(v, cout, NULL,t,G, h);  //Segmentation fault pour penser en cout[V] = +infit*
            double c = ((*U)->cost )+ weight[G.texture[x][y]];

            if(!nodeInHeap(Q,V) ){
              printf("inserting: (%d,%d) dans Q \n",V->pos.x,V->pos.y);
              heap_add(Q,&V);
              G.mark[x][y] = MK_FRONT;
              
            }else if (c >= V->cost){
              continue;
            }
            
            int a = getPos(Q,V);
            printf("a: %d\n", a);
            (*(node *) Q->array[a])->cost = c;
            (*(node *) Q->array[a])->parent = *U;
            (*(node *) Q->array[a])->score = V->cost + h(V->pos,t,&G);

            printf("new score: %f\n",(*(node *) Q->array[a])->score = V->cost + h(V->pos,t,&G) );
        }
      } 
    } 
    compteur+=1;
    
  }
  
  printf("MPT: %d\n", heap_empty(Q));
  /*Troisème Partie*/

  return -1;
}


// Améliorations à faire seulement quand vous aurez bien avancé:
//
// (1) Lorsqu'il y a peu d'obstacles, l'ensemble des sommets visités
//     peut être relativement chaotique (on dirait qu'A* oublie de
//     visiter des sommets), ce qui en fait peut être tout à fait
//     normal. Modifiez la fonction de comparaison pour qu'en cas
//     d'égalité des scores, elle tienne compte des coûts. Observez la
//     différence de comportement si vous privilégiez les coûts
//     croissant ou décroissant.
//
// (2) Le chemin a tendance à zigzaguer, c'est-à-dire à utiliser aussi
//     bien des arêtes horizontales que diagonales (qui peuvent avoir
//     le même coût), même pour des chemins en ligne droite. Essayez
//     de rectifier ce problème d'esthétique en modifiant le calcul de
//     score[v] de sorte qu'à coût[v] égale les arêtes (u,v)
//     horizontales ou verticales soient favorisées (un score plus
//     faible). Il est aussi possible d'intervenir au niveau de la
//     fonction de comparaison des noeuds pour tenir compte de la
//     position du parent. Bien sûr, votre modification ne doit en
//     rien changer la distance (la somme des coûts) entre .start et
//     .end.
//
// (3) Modifier votre implémentation du tas dans heap.c de façon à
//     utiliser un tableau de taille variable, en utilisant realloc()
//     et une stratégie "doublante": lorsqu'il n'y a pas plus assez de
//     place dans le tableau, on double sa taille avec un realloc().
//     On peut imaginer que l'ancien paramètre 'nmax' devienne non pas
//     le nombre maximal d'éléments, mais sa taille maximale initiale
//     (comme par exemple nmax=4).
//
// (4) Gérer plus efficacement la mémoire en libérant les noeuds
//     devenus inutiles. Pour cela on ajoute un champs .nchild à la
//     structure node, permettant de gérer le nombre de fils qu'un
//     node de P ou Q possède. C'est relativement léger et facile à
//     gérer puisqu'on augmente .nchild de u chaque fois qu'on fait
//     parent[v]=p, soit juste après "node v = createNode(p,...)".
//     Pensez à faire .nchild=0 dans createNode(). Notez bien qu'ici
//     on parle de "node", donc de copie de sommet.
//
//     L'observation utile est que tous les nodes de Q sont des
//     feuilles. On va alors pouvoir se débarrasser de tous les nodes
//     ancêtres de ces feuilles simplement en extrayant les nodes de Q
//     dans n'importe quel ordre, par exemple avec une boucle (*) sur
//     heap_pop(Q). On supprime alors chaque node, en mettant à jour
//     le nombre de fils de son père, puis en supprimant le père s'il
//     devient feuille à son tour (son .nchild passant à 0) et ainsi
//     de suite récursivement. On élimine ainsi totalement l'arbre par
//     branches qui se terminent dans Q.
//
//     (*) Si on veut être plus efficace que O(|Q|xlog|Q|), on peut en
//         temps O(|Q|) vider le tableau .array[] directement sans
//         passer par heap_pop(). Mais pour être propre, c'est-à-dire
//         n'utiliser que des primitives de la structure de données et
//         pas son implémentation en accédant à h->array, il est bien
//         vu d'ajouter une primitive "void* heap_tail(heap h)" qui
//         permettrait d'extraire en temps constant le dernier élément
//         du tas (s'il existe). Cela a comme avantage de garder à
//         tout moment la structure de tas.
//
//     Ce processus peut ne pas éliminer tous les nodes car il peut
//     rester des branches de l'arbre qui se terminent par une feuille
//     qui est dans P (et donc qui ont pas de fils dans Q).
//     L'observation est que de tels nodes ne peuvent pas être sur le
//     chemin s->t. On peut donc les supprimer au fur et à mesure
//     directement dans la boucle principale quand on les détecte. On
//     voit qu'un tel node apparaît si après avoir parcouru tous les
//     voisins de u aucun node v n'est créé (et ajouté dans Q). Il
//     suffit donc de savoir si on est passé par heap_add() (ou encore
//     de comparer la taille de Q avant et après la boucle sur les
//     voisins). Si u est une feuille, on peut alors supprimer le node
//     u, mettre à jour .nchild de son père et remonter la branche
//     jusqu'à trouver un node qui n'est plus une feuille. C'est donc
//     la même procédure d'élagage que précédemment qu'on pourrait
//     capturer par une fonction freeNode(node p).


int main(int argc, char *argv[]){

  alpha = (argc>=2)? atof(argv[1]) : 0; // alpha = 0 => Dijkstra par défaut
  unsigned seed=time(NULL)%1000;
  printf("seed: %u\n",seed); // pour rejouer la même grille au cas où
  srandom(seed);

  // testez différentes grilles ...

  grid G = initGridPoints(80,60,TX_FREE,1); // petite grille vide, sans mur
  //grid G = initGridPoints(width,height,TX_FREE,1); // grande grille vide, sans mur
  //grid G = initGridPoints(32,24,TX_WALL,0.2); // petite grille avec quelques murs
  //grid G = initGridLaby(12,9,8); // petit labyrinthe aléatoire
  //grid G = initGridLaby(width/8,height/8,3); // grand labyrinthe aléatoire
  //grid G = initGridFile("mygrid.txt"); // grille à partir d'un fichier modifiable

  // ajoutez à G une (ou plus) "région" de texture donnée ...
  // (déconseillé pour initGridLaby() et initGridFile())

  //addRandomBlob(G, TX_WALL,   (G.X+G.Y)/20);
  //addRandomBlob(G, TX_SAND,   (G.X+G.Y)/15);
  //addRandomBlob(G, TX_WATER,  (G.X+G.Y)/6);
  //addRandomBlob(G, TX_MUD,    (G.X+G.Y)/3);
  //addRandomBlob(G, TX_GRASS,  (G.X+G.Y)/15);
  //addRandomBlob(G, TX_TUNNEL, (G.X+G.Y)/4);
  //addRandomArc(G, TX_WALL,    (G.X+G.Y)/25);

  // sélectionnez des positions s->t ...
  // (inutile pour initGridLaby() et initGridFile())

  G.start=(position){0.1*G.X,0.2*G.Y}, G.end=(position){0.8*G.X,0.9*G.Y};
  //G.start=randomPosition(G,TX_FREE); G.end=randomPosition(G,TX_FREE);

  // constantes à initialiser avant init_SDL_OpenGL()
  scale = fmin((double)width/G.X,(double)height/G.Y); // zoom courant
  hover = false; // interdire les déplacements de points
  init_SDL_OpenGL(); // à mettre avant le 1er "draw"
  drawGrid(G); // dessin de la grille avant l'algo
  update = false; // accélère les dessins répétitifs

  double d = A_star(G,halpha); // heuristique: h0, hvo, alpha×hvo

  // chemin trouvé ou pas ?
  if (d < 0){
    printf("path not found!\n");
    if(G.end.x<0 || G.end.y<0 || G.end.x>=G.X || G.end.y>=G.Y )
      printf("(destination out of the grid)\n");
    if(G.texture[G.end.x][G.end.y] == TX_WALL)
      printf("(destination on a wall)\n");
  }
  else printf("bingo!!! cost of the path: %g\n", d);

  // compte le nombre de sommets explorés pour comparer les
  // heuristiques
  int m = 0;
  for (int i=0; i<G.X; i++)
    for (int j=0; j<G.Y; j++)
      m += (G.mark[i][j] != MK_NULL);
  printf("#nodes explored: %i\n", m);

  while (running) {    // affiche le résultat et attend
    update = true;     // force l'affichage de chaque dessin
    drawGrid(G);       // dessine la grille (/!\ passe update à false)
    handleEvent(true); // attend un évènement
  }

  freeGrid(G);
  cleaning_SDL_OpenGL();
  return 0;
}
