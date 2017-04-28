#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <string.h>
#include <fstream>



#include <vMatrix.h>

using namespace std;


// =============================================================================

bool
retrieveFileDimensions(char* filename,  size_t max_line_size, size_t &nb_lines, size_t &nb_cols);

double**
readMatrixFromFile(char* filename, size_t max_line_size, size_t nb_lines, size_t nb_cols, size_t offset_lines, size_t offset_cols);


// =============================================================================


//fonction acceptable par map
double 
mafonction(double x)
{
  return exp(x + 2) * 3.3f;
}

  

int
main(void)
{
  vMatrix<double> J(4, 3, 1);
  
  cout << "toto" << endl << J;
  J.scalar_mult(2.5f);
  cout << J;

  J.map(&mafonction);
  cout << J;

  double v12 = J[1][2];
  cout << endl << endl << v12;

  J[2][2] = 8.88;
  cout << endl << J;

  vMatrix<double> N(J); 
  cout << N;

  vMatrix<double> O(4, 3, 0.0f);
  O = N;
  cout << endl << O;

  double** mm;
  mm = readMatrixFromFile("mm.txt", 1000, 3, 4, 0, 0);
  vMatrix<double> TM(3, 4, mm);
  cout << endl << TM << endl;
  
  vMatrix<double> tTM(4, 3);
  tTM.transpose(TM);
  cout << endl << tTM << endl;


  vMatrix<double> MMA(3, 2);
  MMA[0][0] = 1;
  MMA[0][1] = 2;
  MMA[1][0] = 3;
  MMA[1][1] = 4;
  MMA[2][0] = 5;
  MMA[2][1] = 6;

  vMatrix<double> MMB(2, 2);
  MMB[0][0] = 5;
  MMB[0][1] = 6;
  MMB[1][0] = 7;
  MMB[1][1] = 8;

  vMatrix<double> MMR(3, 2);
  MMR.multiply(MMA, MMB);
  
  double** t77 = readMatrixFromFile("7x7.txt", 1000, 7, 7, 0, 0);
  vMatrix<double> m77(7, 7, t77);

  double** t27 = readMatrixFromFile("7x7.txt", 1000, 2, 7, 0, 0);
  vMatrix<double> m27(2, 7, t27);

  vMatrix<double> mr(2, 7);
  mr.multiply(m27, m77);

  double** aza = readMatrixFromFile("123.txt", 1000, 3, 3, 0, 0);
  double** azb = readMatrixFromFile("123diag.txt", 1000, 3, 3, 0, 0);
  vMatrix<double> za(3, 3, aza);
  vMatrix<double> zb(3, 3, azb);
  vMatrix<double> zc(3, 3);
  vMatrix<double> zd(3, 3);
  zc.multiply(zb, za);
  cout << endl << zc;

  vMatrix<double> vec(1, 3);
  vec[0][0] = 1; vec[0][1] = 2; vec[0][2] = 3; 
  zd.multiply_from_diag(vec, za);
  cout << endl << zd;
  

  cout << endl << endl << "--------------------------------------------------" << endl;
  cout << m77 << endl << m27 << endl << mr << endl;


  cout << MMA << endl << MMB << endl <<MMR;
  cout << endl << endl << "--------------------------------------------------" << endl;

  // *****************************************************************************


  size_t Wnlig, Wncol;
  size_t Mnlig, Mncol;
  
  
  if (!retrieveFileDimensions("W.txt", 1000, Wnlig, Wncol))
    {cerr << endl << "fichier non carré"; return 1;}

  cout << endl << Wnlig << " " << Wncol << " ";
    
  if (!retrieveFileDimensions("megamat_headless.txt", 1000, Mnlig, Mncol))
    {cerr << endl << "fichier non carré"; return 1;}
  
  cout << endl << Mnlig << " " << Mncol << " ";


  double** w;
  w = readMatrixFromFile("W2.txt", 1000, 7, 7, 1, 1);
  vMatrix<double> W(7, 7, w);
  
  
  // pour voir 16 décimales
  cout.precision(16);
  cout << endl << "W: " << W;

  
  double** megamat;
  Mncol--; Mnlig--;// we don't want 1st line and 1st tokens
  megamat = readMatrixFromFile("megamatc.txt", 1000, Mnlig, Mncol, 1, 1); 
  // 1: skip first token of each line
  vMatrix<double> M(Mnlig, Mncol, megamat);
    

  double** theta;
  theta = readMatrixFromFile("theta2.txt", 1000, 7, 1, 1, 1);
  vMatrix<double> Theta(7, 1, theta);

  //vMatrix<double> tTheta(7, 1);
  //tTheta.transpose(Theta);

  cout << endl << "Theta: " << Theta;
  
  // --------------------------------------------------
  // TA = mA * mB * M
  // TA: {1,7}
  // mA: {1, 30k} que des 1
  // mB: {30k, 30k}
  // mC: {30k, 7}
  
  // --------------------------------------------------
  // M déjà prête, de taille {Mnlig Mncol}
  vMatrix<double> mA(1, Mnlig, 1);
  // mB: trop grosse et creuse: on ne l'alloue pas

  // --------------------------------------------------
  //mB = diag(exp(M * Theta)
  //M * tTheta: {30k, 1}
  vMatrix<double> B(Mnlig, 1);
  vMatrix<double> tB(1, Mnlig);
  B.multiply(M, Theta);
  // visu provisoire
  tB.transpose(B);
  tB.display_limited = true;
  cout << endl << "megamat * theta (t) " << tB << endl;

  B.map(&exp);
  tB.transpose(B);
  

  // --------------------------------------------------
  //mB * M calculée directement => mBM
  vMatrix<double> mBM(Mnlig, 7);
  mBM.multiply_from_diag(tB, M);
  
  // --------------------------------------------------
  // TA = mA * mBM
  vMatrix<double> ma(1, Mnlig, 1);
  vMatrix<double> TA(1, 7);
  TA.multiply(ma, mBM);
  
  // --------------------------------------------------
  // varmu = TA * W * tTA  // {1, 7} * {7, 7} * {7, 1} // scalaire
  vMatrix<double> TAW(1, 7);
  TAW.multiply(TA, W);
  vMatrix<double> tTA(7, 1);
  tTA.transpose(TA);
  vMatrix<double> vM(1, 1);
  vM.multiply(TAW, tTA);
  
  double varmu = vM[0][0];
  double mu = 1849437;
  double ICl = mu - 2*sqrt(varmu);
  double ICh = mu + 2*sqrt(varmu);
  

  // --------------------------------------------------
  // --------------------------------------------------

  cout << endl <<  "--------------------------------------------------" << endl <<  "--------------------------------------------------";

  tB.display_limited = true;
  cout.precision(22);
  //cout << endl << tB << endl;

  M.display_limited = true;
  //cout << M;

  cout << endl << "t1: " << TA;
  cout << endl << "varmu: " << varmu;
  cout << endl << "IC: [" << ICl << " ; " << ICh << " ]" << endl;

  return 0;
}



// =============================================================================


// analyse le fichier texte donné en entrée et détermine le nombre de lignes et de colonnes
// si le nombre de colonnes varie d'une ligne à l'autre, retourne faux
// sinon, écrit résultat dans nb_lines et nb_cols
// indiquer dans max_line_size un nombre de caractères suffisant pour contenir une ligne (typiquement 80, mettre 1000 par précaution)


bool
retrieveFileDimensions(char* filename, size_t max_line_size, size_t &nb_lines, size_t &nb_cols)
{
  char* read_line;
  ARRAY_ALLOC(read_line, max_line_size, char);
  
  ifstream infile (filename);
  assert (infile.is_open());
    
  nb_lines = 0;
  nb_cols = 0;

  char* token;
  size_t nb_tokens = 0;

  while (! infile.eof() )
    {
      infile.getline(read_line, max_line_size, '\n');
      nb_lines++;

      // guess how many tokens are present in line
      nb_tokens = 0;

      token = strtok(read_line, " ");
      while(token != NULL)
	{
	  nb_tokens++;
	  //printf("%s\n", token);
	  token = strtok(NULL, " ");
	}
      
      // assuming all lines have same number of tokens
      if (nb_lines <= 1)
	{nb_cols = nb_tokens;}
      else
	{
	  if (nb_tokens != nb_cols) 
	    if (nb_tokens == 0)
	      { nb_lines--; return true;} // fin de fichier
	    else
	      return false;
	}
      

    }
  
  return true;
}



// =============================================================================


// alloue une matrice et lui affecte les données du fichier spécifié
// matrice de taille nb_lines * nb_cols
// si le fichier contient plus de lignes ou colonnes, elles sont ignorées
// la lecture saute offset_lines linges au début et offset_cols colonnes si besoin (commentaires en première ligne par exemple) 

double**
readMatrixFromFile(char* filename, size_t max_line_size, size_t nb_lines, size_t nb_cols, size_t offset_lines, size_t offset_cols)
{
  char* read_line;
  double** matrix;
  ARRAY_ALLOC(read_line, max_line_size, char);
  MATRIX_ALLOC(matrix, nb_lines, nb_cols, double);

  ifstream infile (filename);
  assert (infile.is_open());

  char* token;
  size_t nb_tokens_read = 0;
  size_t nb_lines_read = 0;
  size_t i = 0, j = 0;


  while (! infile.eof()) 
    {
      infile.getline(read_line, max_line_size, '\n');
      nb_lines_read++;

      if (nb_lines_read <= offset_lines) continue;
      if (nb_lines_read - offset_lines > nb_lines) break;
      
      // guess how many tokens are present in line
      j = 0;

      token = strtok(read_line, " ");
      nb_tokens_read = 1; 
      while(token != NULL)
	{
	  if (nb_tokens_read - offset_cols > nb_cols) break;
	  if (nb_tokens_read > offset_cols) 
	    {
	      //printf("%s\n", token);
	      matrix[i][j] = atof(token);
	      j++;
	    }
	  token = strtok(NULL, " ");
	  nb_tokens_read++;
	}
      
    
      i++;
    }
  
  return matrix;
}
