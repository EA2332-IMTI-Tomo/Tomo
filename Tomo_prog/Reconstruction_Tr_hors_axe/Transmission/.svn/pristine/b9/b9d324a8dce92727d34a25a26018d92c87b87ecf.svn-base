#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "macros.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include <octave-2.9.9/octave/oct.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
//#include "/home/mat/Mat/Projet_C/ImageMagick-6.3.2/Magick++/lib/Magick++.h"
/// #include <Magick++.h>
//#include "/usr/include/ImageMagick/Magick++.h"
#include <fftw3.h>
#include <cstring>
#include <fstream>
#include <sstream>
//#include <strstream>//obsolète!
#include <malloc.h>
//#include "stdafx.h"
#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cv.h>
#include <highgui.h>

typedef struct {
  double Re,Im;
}nb_complexe;
typedef struct {
  int x,y;
}Var2D;
typedef struct {
  int x,y,z;
}Var3D;
typedef struct {
  double Re,Im;
}nbCplx;

using namespace cv;
using namespace std;
///using namespace Magick;

typedef enum PRECISION {DOUBLE, FLOAT, INT, UINT, CHAR};
//pixel information
struct PIXEL
{
  //int x;					//x coordinate of the pixel
  //int y;					//y coordinate
  int increment;			//No. of 2*pi to add to the pixel to unwrap it
  int number_of_pixels_in_group;	//No. of pixels in the pixel group
  float value;			//value of the pixel
  float reliability;
  int group;				//group No.
  int new_group;
  struct PIXEL *head;		//pointer to the first pixel in the group in the linked list
  struct PIXEL *last;		//pointer to the last pixel in the group
  struct PIXEL *next;		//pointer to the next pixel in the group
};


//the EDGE is the line that connects two pixels.
//if we have S PIXELs, then we have S horizental edges and S vertical edges
struct EDGE
{
  float reliab;			//reliabilty of the edge and it depends on the two pixels
  PIXEL *pointer_1;		//pointer to the first pixel
  PIXEL *pointer_2;		//pointer to the second pixel
  int increment;			//No. of 2*pi to add to one of the pixels to unwrap it with respect to the second
};

static float PI = 3.141592654;
static float TWOPI = 6.283185307;
int retroPropag(double *spectre3D_Re,double*spectre3D_Im, double * sup_redon, int dim_final, double *Spectre2D_Re,double *Spectre2D_Im,
                Var2D posSpec, Var3D decal,Var2D NMAX, double rayon);
void changeDim2D(double* tab, double* tabFinal, Var2D dimInit,Var2D dimFin);
void rempli_tableau(unsigned char *finalArray, string path, Var2D coin, Var2D taille);
unsigned char* charger_image2D(unsigned char* phasei,int numero, char* chemin,int coin_x,int coin_y,int taille_x,int taille_y);
void circshift(double* entree,double* entree_shift, int dimx,int decal_x);
void circshift2(double* entree, double* result, Var2D dim,Var2D decal);
void TF2D_1P(double entree[],double sortie[],int taille_x,int taille_y);
void TF2D(double entree_reelle[],double entree_imag[],double fft_reel[],double fft_imag[],int taille_x,int taille_y);
void TF2D_INV(double entree_reelle[],double entree_imag[], double sortie_reelle[],double sortie_imag[],int taille_x,int taille_y);
void TF2D_INV_1P(double entree[],double sortie[],int taille_x,int taille_y);
void tukey2D(double* tuk2D, int dimx,int dimy, float alpha);
void concatener(char chaine1[],char chaine2[],char resultat[]);
void gaussienne(double *tab, int Tx, int sigma, float A, int Exy);
void antigaussienne(double *tab, int Tx, int sigma, float A, int Exy);
void circshift3D(double *volume3D, double *volume3D_shift,int taille_x,int taille_y,int taille_z);
void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D);
void genereCache(double masque[], int t_image, int t_mask, int centreX, int centreY);
void interp3D(double *volume_interp_3D, int taille_x,int taille_y,int taille_z);
void ecrire_rapport(int NXMAX,float rayon,float Rf, float K, int DIMX_CCD2,int coin_x, int coin_y,short int precision_exportation,char *chemin,int nb_proj,float n1,float NA,float Tp, int G);
//void fftw_plan_with_nthreads(int nthreads);
void multiplier_masque(double image[], unsigned char masque[], int t_image, int t_mask, int centreX, int centreY);
void multiplier_masque2(double image[], double masque[], int t_image, int t_mask, int centreX, int centreY);
void crop_window(double* img_src, size_t src_dimx, size_t src_dimy, double* img_dst, size_t dst_dimx, size_t dst_dimy, size_t ul_corner_x, size_t ul_corner_y);
void phase2pi(double* phaseMod1pi, double* partie_reel, double* partie_imag, int taille_x, int taille_y);
void SAV(double *var_sav, int taille, const char *chemin, enum PRECISION precision, char options[]);
void decalCoupe(double *fft_reel, double *fft_imag, double *fft_reel_tmp, double *fft_imag_tmp, Var2D NMAX,Var2D dimCCD);
void phaseUnwrapping(double* partie_reel, double* partie_imag, double* WrappedImage,
                     double* UnwrappedImage, PIXEL *pixel, EDGE *edge, int NXMAX);
void corrplan_I(double *plan_reel, double *plan_imag, double *plan_reel_fit, double *plan_imag_fit, int taille,  Mat masqueC);
void corrPhase(double *phase, double *phase_fit, int taille,  Mat masqueC);
void corrAmpli(double *ampli, double *ampli_fit, int taille,  Mat masqueC);

double Poly2DEval(Mat Coefficients, int deg,int x, int y);
int TaillePoly2D(int deg);
int CountM(Mat& Masque);
void CalculPoly(Mat& imagebrute,Mat& Masque, Mat& polynome,int deg, bool method);
void CalculFond(Mat coefficients, int deg, Mat imagefond);
void read_dataM(char *inputfile,Mat* Data, int image_size, int IMAGE_NUM);
void write_dataM(char *outputfile,Mat Data);
void divComplx(Mat a, Mat b, Mat c, Mat d, Mat& x, Mat& y);
void decal2DGen(double* entree, double* result, Var2D dim,Var2D decal);
Var2D corr_crois(double* obj2D_A, double* obj2D_B, Var2D dim);
void TF2Dcplx(nbCplx *entree, nbCplx *fft, Var2D dim);
void TF2Dcplx_INV(nbCplx *fft_entree, nbCplx *objet, Var2D dim);
void Resharp(double* entree, double* sortie, int taille);
void recalAmpPha(double* amplitude, double* phase, double * ampliRef, double * phaseRef, Var2D dimFin, double* sortie_reel, double* sortie_imag,
                 int cpt_angle, Var2D posDecalAmpli, Var2D posDecalPhase, int NbAngle);
void calculPlan(unsigned char* holo1, unsigned char* holo2, unsigned char* holo3, unsigned char* holo4, double* plan_reel, double* plan_imag,
                int taille, Var2D coin, char* Fichier_holo, double* masque, unsigned char* mod_ref);
void calculPlanCritLG(unsigned char* holo1, unsigned char* holo2, unsigned char* holo3, unsigned char* holo4, double* plan_reel, double* plan_imag,
		      int taille, Var2D coin, char* Fichier_holo, double* masque, unsigned char* mod_ref, int cpt_angle);
void calculPlan1(unsigned char* holo1, unsigned char* holo2, unsigned char* holo3, unsigned char* holo4, double* plan_reel, double* plan_imag,
		 int taille, Var2D coin, char* Fichier_holo, double* masque, unsigned char* mod_ref, double * supRedon, double * supRedon_shift, int dimfinal);

void  retroPropagSA(float deltaZ, double *fft_shift_normRe, double *fft_shift_normIm, double * spectrePropag_normRe, double * spectrePropag_normIm,
                    Var2D dimSpctHolo, double rayon, Var2D CoordSpecI);
void Som_Plan_ds_Vol(double *Vol3DRe, double *Vol3DIm, double *plan2DRe, double *plan2DIm, Var3D dimVol, Var2D dimPlan, int z3Di);
void Grad(double* entreeRe, double* entreeIm, double* Gradient, int taille, int NbPlan, int cptPlan);
void Plan_ds_Vol(double *Vol3DRe, double *Vol3DIm, double *plan2DRe, double *plan2DIm, Var3D dimVol, Var2D dimPlan, int z3Di);

/* conneries 
   void AutofocusGrad(double *fft_shift_normRe, double *fft_shift_normIm, double *sortieRe_refocal, double *sortieIm_refocal,
   Var2D dimFin, Var3D dimVol, double rayon, Var2D CoordSpecI, int cpt_angle);
*/

void AutofocusMod(double *fft_shift_normRe, double *fft_shift_normIm, double *sortieRe_refocal, double *sortieIm_refocal,
		  Var2D dimFin, Var3D dimVol, double rayon, Var2D CoordSpecI, int cpt_angle);
double* alloPlan(int dimx, int dimy);
double* specSynMax(double *planSynRe, double *planSynIm, int NbPlanZ, Var2D dimFin);
void ContrSpeckle(double* holo1, int taille, double contrSpeckle);


///#################################################algo déroulement de phase########################################################################


void AutofocusGrad(double *fft_shift_normRe, double *fft_shift_normIm, double *sortieRe, double *sortieIm, Var2D dimInit, double rayon, Var2D CoordSpecI, int cpt_angle);


void AutofocusGrad_jojo(double *fft_shift_normRe, double *fft_shift_normIm, double *sortieRe, double *sortieIm, size_t dimX, size_t dimY, double rayon, size_t SpecPosX, size_t SpecPosY, int cpt_angle)
{
  
  Var2D size = {(int)dimX, (int)dimY};
  Var2D specPos = {(int)SpecPosX, (int)SpecPosY};
  AutofocusGrad(fft_shift_normRe, fft_shift_normIm, sortieRe, sortieIm, size, rayon, specPos, cpt_angle);
  
}


///#################################################algo déroulement de phase########################################################################


//another version of Mixtogether but this function should only be use with the sort program
void  Mix(EDGE *Pointer1, int *index1, int *index2, int size)
{
  int counter1 = 0;
  int counter2 = 0;
  int *TemporalPointer = index1;

  int *Result = (int *) calloc(size * 2, sizeof(int));
  int *Follower = Result;

  while ((counter1 < size) && (counter2 < size))
    {
      if ((Pointer1[*(index1 + counter1)].reliab <= Pointer1[*(index2 + counter2)].reliab))
	{
	  *Follower = *(index1 + counter1);
	  Follower++;
	  counter1++;
	}
      else
        {
	  *Follower = *(index2 + counter2);
	  Follower++;
	  counter2++;
        }
    }//while

  if (counter1 == size)
    {
      memcpy(Follower, (index2 + counter2), sizeof(int)*(size-counter2));
    }
  else
    {
      memcpy(Follower, (index1 + counter1), sizeof(int)*(size-counter1));
    }

  Follower = Result;
  index1 = TemporalPointer;

  int i;
  for (i=0; i < 2 * size; i++)
    {
      *index1 = *Follower;
      index1++;
      Follower++;
    }

  free(Result);
}

///-----------------------start quick_sort algorithm --------------------------------

//this is may be the fastest sort program;
//see the explination in quickSort function below
void  sort(EDGE *Pointer, int *index, int size)
{
  if (size == 2)
    {
      if ((Pointer[*index].reliab) > (Pointer[*(index+1)].reliab))
	{
	  int Temp;
	  Temp = *index;
	  *index = *(index+1);
	  *(index+1) = Temp;
	}
    }
  else if (size > 2)
    {
      sort(Pointer, index, size/2);
      sort(Pointer, (index + (size/2)), size/2);
      Mix(Pointer, index, (index + (size/2)), size/2);
    }
}

//this function tries to implement a nice idea explained below
//we need to sort edge array. Each edge element conisists of 16 bytes.
//In normal sort program we compare two elements in the array and exchange
//their place under some conditions to do the sorting. It is very probable
// that an edge element may change its place hundred of times which makes
//the sorting a very time consuming operation. The idea in this function
//is to give each edge element an index and move the index not the edge
//element. The edge need 4 bytes which makes the sorting operation faster.
// After finishingthe sorting of the indexes, we know the position of each index.
//So we know how to sort edges
void  quick_sort(EDGE *Pointer, int size)
{
  int *index = (int *) calloc(size, sizeof(int));
  int i;

  for (i=0; i<size; ++i)
    index[i] = i;

  sort(Pointer, index, size);

  EDGE * a = (EDGE *) calloc(size, sizeof(EDGE));
  for (i=0; i<size; ++i)
    a[i] = Pointer[*(index + i)];

  memcpy(Pointer, a, size*sizeof(EDGE));

  free(index);
  free(a);
}

///-------------------------end quick_sort algorithm --------------------------------


///---------------------start quicker_sort algorithm --------------------------------
#define swap(x,y) {EDGE t; t=x; x=y; y=t;}
#define order(x,y) if (x.reliab > y.reliab) swap(x,y)
#define o2(x,y) order(x,y)
#define o3(x,y,z) o2(x,y); o2(x,z); o2(y,z)

typedef enum {yes, no} yes_no;

yes_no find_pivot(EDGE *left, EDGE *right, float *pivot_ptr)
{
  EDGE a, b, c, *p;

  a = *left;
  b = *(left + (right - left) /2 );
  c = *right;
  o3(a,b,c);

  if (a.reliab < b.reliab)
    {
      *pivot_ptr = b.reliab;
      return yes;
    }

  if (b.reliab < c.reliab)
    {
      *pivot_ptr = c.reliab;
      return yes;
    }

  for (p = left + 1; p <= right; ++p)
    {
      if (p->reliab != left->reliab)
	{
	  *pivot_ptr = (p->reliab < left->reliab) ? left->reliab : p->reliab;
	  return yes;
	}
      return no;
    }
}

EDGE *partition(EDGE *left, EDGE *right, float pivot)
{
  while (left <= right)
    {
      while (left->reliab < pivot)
	++left;
      while (right->reliab >= pivot)
	--right;
      if (left < right)
	{
	  swap (*left, *right);
	  ++left;
	  --right;
	}
    }
  return left;
}

void quicker_sort(EDGE *left, EDGE *right)
{
  EDGE *p;
  float pivot;

  if (find_pivot(left, right, &pivot) == yes)
    {
      p = partition(left, right, pivot);
      quicker_sort(left, p - 1);
      quicker_sort(p, right);
    }
}

///--------------end quicker_sort algorithm -----------------------------------

//--------------------start initialse pixels ----------------------------------
//initialse pixels. See the explination of the pixel class above.
//initially every pixel is a gorup by its self
void  initialisePIXELs(double *WrappedImage, PIXEL *pixel, int image_width, int image_height)
{
  PIXEL *pixel_pointer = pixel;
  double *wrapped_image_pointer = WrappedImage;
  int i, j;

  for (i=0; i < image_height; i++)
    {
      for (j=0; j < image_width; j++)
        {
	  //pixel_pointer->x = j;
	  //pixel_pointer->y = i;
	  pixel_pointer->increment = 0;
	  pixel_pointer->number_of_pixels_in_group = 1;
	  pixel_pointer->value = *wrapped_image_pointer;
	  pixel_pointer->reliability = 9999999+rand();
	  pixel_pointer->head = pixel_pointer;
	  pixel_pointer->last = pixel_pointer;
	  pixel_pointer->next = NULL;
	  pixel_pointer->new_group = 0;
	  pixel_pointer->group = -1;
	  pixel_pointer++;
	  wrapped_image_pointer++;
	}
    }
}
//-------------------end initialise pixels -----------

//gamma function in the paper
float wrap(float pixel_value)
{
  float wrapped_pixel_value;
  if (pixel_value > PI)	wrapped_pixel_value = pixel_value - TWOPI;
  else if (pixel_value < -PI)	wrapped_pixel_value = pixel_value + TWOPI;
  else wrapped_pixel_value = pixel_value;
  return wrapped_pixel_value;
}

// pixelL_value is the left pixel,	pixelR_value is the right pixel
int find_wrap(float pixelL_value, float pixelR_value)
{
  float difference;
  int wrap_value;
  difference = pixelL_value - pixelR_value;

  if (difference > PI)	wrap_value = -1;
  else if (difference < -PI)	wrap_value = 1;
  else wrap_value = 0;

  return wrap_value;
}

void calculate_reliability(double *wrappedImage, PIXEL *pixel, int image_width, int image_height)
{
  int image_width_plus_one = image_width + 1;
  int image_width_minus_one = image_width - 1;
  PIXEL *pixel_pointer = pixel + image_width_plus_one;
  double *WIP = wrappedImage + image_width_plus_one; //WIP is the wrapped image pointer
  float H, V, D1, D2;
  int i, j;

  for (i = 1; i < image_height -1; ++i)
    {
      for (j = 1; j < image_width - 1; ++j)
	{
	  H = wrap(*(WIP - 1) - *WIP) - wrap(*WIP - *(WIP + 1));
	  V = wrap(*(WIP - image_width) - *WIP) - wrap(*WIP - *(WIP + image_width));
	  D1 = wrap(*(WIP - image_width_plus_one) - *WIP) - wrap(*WIP - *(WIP + image_width_plus_one));
	  D2 = wrap(*(WIP - image_width_minus_one) - *WIP) - wrap(*WIP - *(WIP + image_width_minus_one));
	  pixel_pointer->reliability = H*H + V*V + D1*D1 + D2*D2;
	  pixel_pointer++;
	  WIP++;
	}
      pixel_pointer += 2;
      WIP += 2;
    }
}

//calculate the reliability of the horizental edges of the image
//it is calculated by adding the reliability of pixel and the relibility of
//its right neighbour
//edge is calculated between a pixel and its next neighbour
void  horizentalEDGEs(PIXEL *pixel, EDGE *edge, int image_width, int image_height)
{
  int i, j;
  EDGE *edge_pointer = edge;
  PIXEL *pixel_pointer = pixel;

  for (i = 0; i < image_height; i++)
    {
      for (j = 0; j < image_width - 1; j++)
	{
	  edge_pointer->pointer_1 = pixel_pointer;
	  edge_pointer->pointer_2 = (pixel_pointer+1);
	  edge_pointer->reliab = pixel_pointer->reliability + (pixel_pointer + 1)->reliability;
	  edge_pointer->increment = find_wrap(pixel_pointer->value, (pixel_pointer + 1)->value);
	  pixel_pointer++;
	  edge_pointer++;
	}
      pixel_pointer++;
    }
}

//calculate the reliability of the vertical EDGEs of the image
//it is calculated by adding the reliability of pixel and the relibility of
//its lower neighbour in the image.
void  verticalEDGEs(PIXEL *pixel, EDGE *edge, int image_width, int image_height)
{
  int i, j;

  PIXEL *pixel_pointer = pixel;
  EDGE *edge_pointer = edge + (image_height) * (image_width - 1);

  for (i=0; i<image_height - 1; i++)
    {
      for (j=0; j < image_width; j++)
	{
	  edge_pointer->pointer_1 = pixel_pointer;
	  edge_pointer->pointer_2 = (pixel_pointer + image_width);
	  edge_pointer->reliab = pixel_pointer->reliability + (pixel_pointer + image_width)->reliability;
	  edge_pointer->increment = find_wrap(pixel_pointer->value, (pixel_pointer + image_width)->value);
	  pixel_pointer++;
	  edge_pointer++;
	} //j loop
    } // i loop
}

//gather the pixels of the image into groups
void  gatherPIXELs(EDGE *edge, int image_width, int image_height)
{
  int k;

  //Number of rialiable edges (not at the borders of the image)
  int no_EDGEs = (image_width - 1) * (image_height) + (image_width) * (image_height - 1);
  PIXEL *PIXEL1;
  PIXEL *PIXEL2;

  PIXEL *group1;
  PIXEL *group2;
  EDGE *pointer_edge = edge;
  int incremento;

  for (k = 0; k < no_EDGEs; k++)
    {
      PIXEL1 = pointer_edge->pointer_1;
      PIXEL2 = pointer_edge->pointer_2;

      //PIXEL 1 and PIXEL 2 belong to different groups
      //initially each pixel is a group by it self and one pixel can construct a group
      //no else or else if to this if
      if (PIXEL2->head != PIXEL1->head)
	{
	  //PIXEL 2 is alone in its group
	  //merge this pixel with PIXEL 1 group and find the number of 2 pi to add
	  //to or subtract to unwrap it
	  if ((PIXEL2->next == NULL) && (PIXEL2->head == PIXEL2))
	    {
	      PIXEL1->head->last->next = PIXEL2;
	      PIXEL1->head->last = PIXEL2;
	      (PIXEL1->head->number_of_pixels_in_group)++;
	      PIXEL2->head=PIXEL1->head;
	      PIXEL2->increment = PIXEL1->increment-pointer_edge->increment;
	    }

	  //PIXEL 1 is alone in its group
	  //merge this pixel with PIXEL 2 group and find the number of 2 pi to add
	  //to or subtract to unwrap it
	  else if ((PIXEL1->next == NULL) && (PIXEL1->head == PIXEL1))
	    {
	      PIXEL2->head->last->next = PIXEL1;
	      PIXEL2->head->last = PIXEL1;
	      (PIXEL2->head->number_of_pixels_in_group)++;
	      PIXEL1->head = PIXEL2->head;
	      PIXEL1->increment = PIXEL2->increment+pointer_edge->increment;
	    }

	  //PIXEL 1 and PIXEL 2 both have groups
	  else
            {
	      group1 = PIXEL1->head;
	      group2 = PIXEL2->head;
	      //the no. of pixels in PIXEL 1 group is large than the no. of PIXELs
	      //in PIXEL 2 group.   Merge PIXEL 2 group to PIXEL 1 group
	      //and find the number of wraps between PIXEL 2 group and PIXEL 1 group
	      //to unwrap PIXEL 2 group with respect to PIXEL 1 group.
	      //the no. of wraps will be added to PIXEL 2 grop in the future
	      if (group1->number_of_pixels_in_group > group2->number_of_pixels_in_group)
		{
		  //merge PIXEL 2 with PIXEL 1 group
		  group1->last->next = group2;
		  group1->last = group2->last;
		  group1->number_of_pixels_in_group = group1->number_of_pixels_in_group + group2->number_of_pixels_in_group;
		  incremento = PIXEL1->increment-pointer_edge->increment - PIXEL2->increment;
		  //merge the other pixels in PIXEL 2 group to PIXEL 1 group
		  while (group2 != NULL)
		    {
		      group2->head = group1;
		      group2->increment += incremento;
		      group2 = group2->next;
		    }
		}

	      //the no. of PIXELs in PIXEL 2 group is large than the no. of PIXELs
	      //in PIXEL 1 group.   Merge PIXEL 1 group to PIXEL 2 group
	      //and find the number of wraps between PIXEL 2 group and PIXEL 1 group
	      //to unwrap PIXEL 1 group with respect to PIXEL 2 group.
	      //the no. of wraps will be added to PIXEL 1 grop in the future
	      else
                {
		  //merge PIXEL 1 with PIXEL 2 group
		  group2->last->next = group1;
		  group2->last = group1->last;
		  group2->number_of_pixels_in_group = group2->number_of_pixels_in_group + group1->number_of_pixels_in_group;
		  incremento = PIXEL2->increment + pointer_edge->increment - PIXEL1->increment;
		  //merge the other pixels in PIXEL 2 group to PIXEL 1 group
		  while (group1 != NULL)
		    {
		      group1->head = group2;
		      group1->increment += incremento;
		      group1 = group1->next;
		    } // while

                } // else
            } //else
        } ;//if

      pointer_edge++;
    }
}

//unwrap the image
void  unwrapImage(PIXEL *pixel, int image_width, int image_height)
{
  int i;
  int image_size = image_width * image_height;
  PIXEL *pixel_pointer=pixel;

  for (i = 0; i < image_size; i++)
    {
      pixel_pointer->value += TWOPI * (float)(pixel_pointer->increment);
      pixel_pointer++;
    }
}

//the input to this unwrapper is an array that contains the wrapped phase map.
//copy the image on the buffer passed to this unwrapper to over write the unwrapped
//phase map on the buffer of the wrapped phase map.
void  returnImage(PIXEL *pixel, double *unwrappedImage, int image_width, int image_height)
{
  int i;
  int image_size = image_width * image_height;
  double *unwrappedImage_pointer = unwrappedImage;
  PIXEL *pixel_pointer = pixel;

  for (i=0; i < image_size; i++)
    {
      *unwrappedImage_pointer = pixel_pointer->value;
      pixel_pointer++;
      unwrappedImage_pointer++;
    }
}

void read_data(char *inputfile,double *Data, int length)
{
  printf("Reading the Wrapped Values form Binary File.............>");
  float precision;
  FILE *ifptr;
  ifptr = fopen(inputfile,"rb");
  if(ifptr == NULL) printf("Error opening the file\n");
  for(int cpt=0;cpt<length;cpt++)
    {
      fread(&precision,sizeof(precision),1,ifptr);
      Data[cpt]=precision;
    }
  fclose(ifptr);
  printf(" Done.\n");
}
/*
  void write_data(char *outputfile,double *Data,int length)
  {
  //printf("Writing the Unwrapped Values to Binary File.............>");
  FILE *ifptr;
  ifptr = fopen(outputfile,"a+b");
  if(ifptr == NULL) printf("Error opening the file\n");
  fwrite(Data,sizeof(double),length,ifptr);
  fclose(ifptr);
  //printf(" Done.\n");
  }
*/

void write_data(char *outputfile,double *Data,int length)
{
  //printf("Writing the Unwrapped Values to Binary File.............>");
  FILE *ifptr;
  ifptr = fopen(outputfile,"a+b");
  if(ifptr == NULL) printf("Error opening the file\n");
  for(int cpt=0;cpt<length;cpt++)
    {
      float precision=Data[cpt];
      fwrite(&precision,sizeof(precision),1,ifptr);
    }

  fclose(ifptr);
  //printf(" Done.\n");
}

///#######################################################fonction classique##########################################################################
int cplx_mult(double* Re1,double* Im1,double* Re2,double* Im2,double* resultRe,double* resultIm,Var2D dim)
{
  for(int cpt=0;cpt<dim.x*dim.y;cpt++)
    {
      resultRe[cpt]=Re1[cpt]*Re2[cpt]-Im1[cpt]*Im2[cpt];
      resultIm[cpt]=Re1[cpt]*Im2[cpt]+Im1[cpt]*Re2[cpt];
    }
  return 1;
}

int CreerZoneFresnel(double *FresnelRe,double * FresnelIm, Var2D dim, Var2D centre,float d, float lambda)
{
  float Tp=80*pow(10,-9);
  for(int cpty=0;cpty<dim.y;cpty++)
    {
      int hauteur=cpty*dim.x;
      for(int cptx=0;cptx<dim.x;cptx++)
	{
	  int cpt=hauteur+cptx;
	  FresnelRe[cpt]=cos(PI*Tp*Tp*(pow(cptx-centre.x,2)+pow(cpty-centre.y,2))/(lambda*d));
	  FresnelIm[cpt]=sin(PI*Tp*Tp*(pow(cptx-centre.x,2)+pow(cpty-centre.y,2))/(lambda*d));
	}
    }
  return 1;
}
//Spectre2D_Re=fft_Re_shift
int retroPropag(double *spectre3D_Re,double*spectre3D_Im, double * sup_redon, int dim_final, double *Spectre2D_Re,double *Spectre2D_Im,
                Var2D posSpec, Var3D decal,Var2D NMAX, double rayon)
{
  int k_T=0, k_R=0;
  int ptsExcluHolo=0;
  int dimPlanFinal=round(dim_final*dim_final),
    xm0=posSpec.x,ym0=posSpec.y,
    dimVolX=round(dim_final);
  ///--------création de variable pour éviter N calculs dans la boucle sur le volume 3D
  double r2=rayon*rayon;
  double arg_z_arc=0,z_arc=0,zm0;
  //printf("round(rayon*rayon-(xm0)^2-(ym0)^2: %i\n",round(rayon*rayon-(xm0)^2-(ym0)^2));
  double zm0_carre = rayon*rayon-xm0*xm0-ym0*ym0;
  if(round(zm0_carre)>-1)
    {
      zm0=sqrt(zm0_carre);
      int NXMAX_CARRE=NMAX.x*NMAX.x;
      for (int y = -NMAX.y; y < NMAX.y; y++)
	{ //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
	  int y_carre=y*y;
	  for (int x = -NMAX.x; x < NMAX.x; x++)
	    { //on balaye l'image 2D en y, centre au milieu
	      int cpt=(y+NMAX.y)*2*NMAX.x+x+NMAX.x;//calcul du cpt du tableau 1D de l'image 2D
	      if(x*x+y_carre<NXMAX_CARRE)
		{ //ne pas depasser l'ouverture numérique pour 1 hologramme
		  double z_carre=r2-x*x-y_carre; //altitude au carré des données
		  double z_T=round(sqrt(z_carre)-zm0);
		  double altitude_T=(z_T+decal.z)*dimPlanFinal; //donne n'importequoi sans l'arrondi sur z!!
		  k_T=(-xm0+x+decal.x)+(-ym0+y+decal.y)*dimVolX+round(altitude_T);//indice du tableau 1D du volume 3D
		  spectre3D_Re[k_T]+=Spectre2D_Re[cpt];//pour calculer l'image;//
		  spectre3D_Im[k_T]+=Spectre2D_Im[cpt];//pour calculer l'image
		  sup_redon[k_T]+=1;//pour calculer le support

		  //                                                                double z_R=round(sqrt(z_carre)+zm0);
		  //                                                                double altitude_R=(z_R+decal.z)*dimPlanFinal;
		  //                                                                k_R=(-xm0+x+decal.x)+(-ym0+y+decal.y)*dimVolX+round(altitude_R);
		  //                                                                spectre3D_Re[k_R]+=Spectre2D_Re[cpt];//pour calculer l'image;//
		  //                                                                spectre3D_Im[k_R]+=Spectre2D_Im[cpt];//pour calculer l'image
		  //                                                                sup_redon[k_R]+=1;//pour calculer le support
		}
	      else
		ptsExcluHolo++;
	    }
	}

    }
  return ptsExcluHolo;
}

void decalCoupe(double *fft_reel, double *fft_imag, double *fft_reel_tmp, double *fft_imag_tmp, Var2D NMAX,Var2D dimCCD)
{
  for (int xi=0; xi<NMAX.x; xi++) {
    for (int yi=0; yi<NMAX.y; yi++) {
      int cpt1=yi*dimCCD.x+xi;
      int cpt2=yi*(NMAX.x)*2+xi;

      fft_reel[cpt2]=fft_reel_tmp[cpt1];
      fft_imag[cpt2]=fft_imag_tmp[cpt1];
    }
    for (int yi=dimCCD.y-NMAX.y; yi<dimCCD.y; yi++) {
      int cpt1=yi*dimCCD.x+xi;

      int cpt2 = xi+(yi-dimCCD.y+2*NMAX.y)*2*(NMAX.x);
      fft_reel[cpt2]=fft_reel_tmp[cpt1];
      fft_imag[cpt2]=fft_imag_tmp[cpt1];
    }
  }
  ///---////////////////////////////////////////////deuxieme demi-espace
  for(int xi=dimCCD.x-NMAX.x; xi<dimCCD.x; xi++) {
    for (int yi=0; yi<NMAX.y; yi++) {
      int cpt1=yi*dimCCD.x+xi;
      int cpt2=yi*(NMAX.x)*2+(xi-dimCCD.x+2*NMAX.x);

      fft_reel[cpt2]=fft_reel_tmp[cpt1];
      fft_imag[cpt2]=fft_imag_tmp[cpt1];
    }

    for (int yi=dimCCD.y-NMAX.y; yi<dimCCD.y; yi++) {
      int cpt1=yi*dimCCD.x+xi;

      int cpt2 = (xi-dimCCD.x+2*NMAX.x)+(yi-dimCCD.y+2*NMAX.y)*2*(NMAX.x);
      fft_reel[cpt2]=fft_reel_tmp[cpt1];
      fft_imag[cpt2]=fft_imag_tmp[cpt1];
    }
  }
}
void SAV(double *var_sav, int taille, const char *chemin, enum PRECISION precision, char options[])
{
  FILE *fichier_ID;
  fichier_ID= fopen(chemin, options);
  if(fichier_ID==0)
    cout<<"Erreur d'ouverture du fichier "<<chemin<<endl;

  switch(precision) {
  case DOUBLE: //64 bit

    for(unsigned int cpt=0; cpt<taille; cpt++) {
      double tampon=var_sav[cpt];
      fwrite(&tampon,sizeof(tampon),1,fichier_ID);
    }
    break;
  case FLOAT://32 bits float

    for(unsigned int cpt=0; cpt<taille; cpt++) {
      float tampon=var_sav[cpt];
      fwrite(&tampon,sizeof(tampon),1,fichier_ID);
    }
    break;

  case INT: //32 bit signé

    for(unsigned int cpt=0; cpt<taille; cpt++) {
      int tampon=var_sav[cpt];
      fwrite(&tampon,sizeof(tampon),1,fichier_ID);
    }
    break;
  case UINT://32 bit non signé

    for(unsigned int cpt=0; cpt<taille; cpt++) {
      unsigned int tampon=var_sav[cpt];
      fwrite(&tampon,sizeof(tampon),1,fichier_ID);
    }
    break;
  case CHAR: //8 bits

    for(unsigned int cpt=0; cpt<taille; cpt++) {
      char tampon=var_sav[cpt];
      fwrite(&tampon,sizeof(tampon),1,fichier_ID);
    }
    break;
  default:
    break;
  }

  fclose(fichier_ID);
}
void changeDim2D(double* tab, double* tabFinal, Var2D dimInit, Var2D dimFin)
{
  int diffX=round(dimFin.x-dimInit.x)/2;
  int diffY=round(dimFin.y-dimInit.y)/2;

  for(int x=0;x<dimInit.x;x++)
    {
      for( int y=0;y<dimInit.y;y++)
        {   int cpt_init=y*dimInit.x+x;
	  int cpt_final=(y+diffY)*dimFin.x+x+diffX;
	  tabFinal[cpt_final]=tab[cpt_init];
        }
    }
}


void concatener(char chaine1[], char chaine2[], char resultat[])//fonction limitée mais intéret car copie locale (pas besoin de tampon)
{

  //récupérer la longueur du resultat
  size_t LongrChaine1=strlen(chaine1);
  size_t LongrChaine2=strlen(chaine2);
  short unsigned int LongrChaine=LongrChaine1+LongrChaine2;
  char tampon[LongrChaine];//+1 pour zero final?

  strcpy(tampon,chaine1);
  strcat(tampon,chaine2);

  for(short unsigned int cpt=0; cpt<LongrChaine+1; cpt++)
    resultat[cpt]=tampon[cpt];

}

///#######masque pour &écraser jumeau############""
void antigaussienne(double *tab, int Tx, int sigma, float A, int Exy)
{

  int corr_paire=0;
  if(Tx%2==0)
    {
      cout<<"taille Tx="<<Tx<<" paire : le masque ne serait pas centré!"<<endl;
      corr_paire=1;
    }

  if(sigma==0)
    sigma=1;
  short  int x,y, Tinf=-round(Tx/2),Tsup=round(Tx/2);
  short unsigned int cptx,cpty;
  float Ex,Ey,sigmax,sigmay;
  Ex=Exy;
  Ey=Exy;
  sigmax=sigma;
  sigmay=sigma;
  if(Tx==1)
    tab[0]=0;
  else
    {
      for(x=Tinf;x<Tsup+1-corr_paire;x++)
	{
	  cptx=x+Tsup;
	  for( y=Tinf;y<Tsup+1-corr_paire;y++)
	    {
	      cpty=y+Tsup;
	      tab[cpty*Tx+cptx]=1-A*exp(-(pow((x-Ex),4)/(2*sigmax*sigmax)+pow((y-Ey),4)/(2*sigmay*sigmay)));
	    }
	}
    }

}

void gaussienne(double *tab, int Tx, int sigma, float A, int Exy)
{

  int corr_paire=0;
  if(Tx%2==0)
    {
      cout<<"taille Tx="<<Tx<<" paire : le masque ne serait pas centré!"<<endl;
      corr_paire=1;
    }

  if(sigma==0)
    sigma=1;
  short  int x,y, Tinf=-round(Tx/2),Tsup=round(Tx/2);
  short unsigned int cptx,cpty;
  float Ex,Ey,sigmax,sigmay;
  Ex=Exy;
  Ey=Exy;
  sigmax=sigma;
  sigmay=sigma;
  if(Tx==1)
    tab[0]=0;
  else
    {
      for(x=Tinf;x<Tsup+1-corr_paire;x++)
	{
	  cptx=x+Tsup;
	  for( y=Tinf;y<Tsup+1-corr_paire;y++)
	    {
	      cpty=y+Tsup;
	      tab[cpty*Tx+cptx]=A*exp(-(pow((x-Ex),4)/(2*sigmax*sigmax)+pow((y-Ey),4)/(2*sigmay*sigmay)));
	    }
	}
    }

}

///////////////////////////////////////////////générer un masque tukey/////////////////////////////////////////////////////
void tukey2D(double *tuk2D,int dimx, int dimy, float alpha)

{
  double *tukey1D=new double[dimx];

  int N=dimx;

  const float
    interva=round(alpha*(N-1)/2),
    intervb=round((N-1)*(1-alpha/2));


  for(int i=0;i<interva+1;i++)
    {
      tukey1D[i]=((1+cos(PI*(2*i/(alpha*(N-1))-1)))*0.5);
    }
  for(int i=interva+1;i<intervb+1;i++)
    {
      tukey1D[i]=1;
    }
  for(int i=intervb+1;i<N;i++)
    {
      tukey1D[i]=((1+cos(PI*(2*i/(alpha*(N-1))+1-2/alpha)))*0.5);
    }

  //for (int i=0; i< N; i++)
  //{
  //    cout << tukey1D[i] << endl ;
  //}


  //for (int i=0;i<N*N;i++)
  //{
  //        int ix=i%(N);
  //        int iy=i/(N);
  //        tukey2D[i]= tukey1D[ix] * tukey1D[iy];
  //        cout << tukey2D[i] <<endl;
  //}

  for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
	{
	  int cpt=j*N+i;
	  tuk2D[cpt]= tukey1D[i] * tukey1D[j];
	  //cout << tukey2D[cpt] <<endl;
	}
    }
  delete[] tukey1D;
}

void ecrire_rapport(int NXMAX,float rayon,float Rf, float K, int DIMX_CCD2,int coin_x, int coin_y,short int precision_exportation,char* chemin,int nb_proj,float n1,float NA,float Tp, int G)
{
  time_t date;
  time(&date);
  FILE *fichier_rapport ;
  /*  ouverture pour ecriture (w) en mode texte (t) */
  fichier_rapport = fopen ("/ramdisk/rapport_calcul.txt", "wt") ;
  if (fichier_rapport == NULL)
    {
      printf ("impossible de créer le fichier rapport_calcul.txt\n") ;
      exit (0) ;
    }
  //fprintf(fichier_rapport,"Date     : %s\nNXMAX    : %i\n,Rayon : %i\nRf       : %f\nK : %f\ndimx_ccd : %d\ncoin_x   : %d\ncoin_y   : %d\nPrecision: %i bits\nSession  : %s\nnb_proj  : %i\nindice n1: %f\nNA       : %f\nT_pixel  : %e\nG        : %i",ctime(&date),NXMAX,rayon,Rf,K,DIMX_CCD2,coin_x,coin_y,8*sizeof(precision),chemin,nb_proj,n1,NA,Tp,G);
  fprintf(fichier_rapport,"Date : %s\n NXMAX=%i\n Rf=%f\n Precision: %i bits\n Session  : %s\n Nombre de projections : %i",ctime(&date),NXMAX,Rf,8*precision_exportation,chemin,nb_proj);
  fclose (fichier_rapport);
}


void genereCache(double masque[], int t_image, int t_mask, int centreX, int centreY)
{
  int t_imageX=t_image;
  int t_imageY=t_image;
  if((t_mask%2)!=0)
    cout<<"fonction genere_masque : attention, la dimension de votre masque est impaire! t_mask="<<t_mask<<endl;
  for(int pixel=0;pixel<2*t_image*2*t_image;pixel++)
    {
      masque[pixel]=0;
    }
  for(int x=centreX-t_mask/2;x<centreX+t_mask/2;x++)
    {
      for(int y=centreY-t_mask/2;y<centreY+t_mask/2;y++)
	{   //le masque repasse du coté gauche lorsque le cache touche le bord droit! A corriger (?)
	  int x_2=x;
	  int y_2=y;
	  if(x>2*t_imageX)
	    x_2=x-2*t_imageX;
	  if(y>2*t_imageY)
	    y_2=y-2*t_imageY;
	  if(x<0)
	    x_2=2*t_imageX+x;
	  if(y<0)
	    y_2=2*t_imageY+y;
	  ///coordonnées 1D du centre
	  int cptj=2*t_imageY*y_2+x_2;
	  masque[cptj]=1;
	}
    }
}

/////////////////////////////////////////////////////////////////////////////////////
//fonction chargeant une image 2D. Retourne un pointeur sur le tableau 1D qui la contient
//  Nécessite le nom du pointeur à remplir,
//le numéro du pas de phase shifting (1,2 3 ou 4), le chemin, et cpt_fichier
//pour boucler sur les 1000 projections
/////////////////////////////////////////////////////////////////////////////////////

unsigned char* charger_image2D(unsigned char* phasei,int numero,char *Chemin,int coin_x,int coin_y,int taille_x,int taille_y)
{
  Var2D taille= {taille_x,taille_y}, coin= {coin_x,coin_y};
  char CopieChemin[100];
  char FinNom[15];
  strcpy(CopieChemin,Chemin);//sauvegarder chemin car sprintf écrase
  //sprintf(FinNom,"-00%i.bmp",numero);   // transmission
  sprintf(FinNom,"_p%i.pgm",numero);    //reflexion
  //créer le chemin final
  strcat(CopieChemin, FinNom);
  //printf("Copiechemin vaut : %s\n", CopieChemin);
  //charger l'image dans le tableau
  ///tester la validité du chemin "copiechemin"
  /*FILE *test_chemin = fopen(CopieChemin,"r");

    if (test_chemin == NULL)
    {
    fprintf( stderr, "attention, Ficher %i absent.\n",cpt_fichier); //fichier invalide? alors on retourne une erreur
    phasei=NULL;
    return phasei;
    fclose(test_chemin);
    }
    else //sinon on peut essayer de remplir le tableau
    {
    phasei = rempli_tableau(CopieChemin, coin_x, coin_y,taille_x,taille_y);
    return phasei;
    fclose(test_chemin);
    }*/
  rempli_tableau(phasei, CopieChemin, coin,taille);
  return phasei;
}

/////////////////////////////////////////////////////////////////////////////////////
//Fonction transferant les pixels de l'image 2D de taille (nx,ny) dans un tableau 1D de taille (nx*ny)
// nécessite la bibliotheque imageMagick
/////////////////////////////////////////////////////////////////////////////////////
void rempli_tableau(unsigned char *finalArray, string path, Var2D coin, Var2D taille)
{

  MSG_ASSERT(false, "pas possible");
  /*
    Image Monimage;
    //int i, currentImageWidth, currentImageHeight;

    Monimage.read(path);////// chargement en memoire de l'image
    Monimage.type( GrayscaleType );	////// Mon image est N&B
    //Monimage.display();
    //Monimage.channelDepth();
    Monimage.compressType(NoCompression);
    //Monimage.crop(Geometry(dimx_ccd,dimy_ccd, 0, 0) );
    //Monimage.write("/home/mat/tomo_test/test.bmp");
    //currentImageWidth = Monimage.columns();////// extraction des caractéristiques de l'image
    //currentImageHeight = Monimage.rows();
    //	finalArray = new unsigned char[taille_x*taille_y];////// reservation de la taille memoire du tableau final
    ////// lecture de l'image
    Monimage.getPixels(coin.x,coin.y,taille.x,taille.y);
    ////// ecriture de l'image dans un tableau
    Monimage.writePixels(GrayQuantum, finalArray);
    //printf("valeur de la case 256 %d \n", finalArray[256]);//affichage pour contrôle
    */
  
}

/////////////////////////////////////////////////////////////////////////////////////
//fonction de calcul TF2D
/////////////////////////////////////////////////////////////////////////////////////
void TF2D(double entree_reelle[],double entree_imag[],double fft_reel[],double fft_imag[],int taille_x,int taille_y)
{
  fftw_plan_with_nthreads(6);

  int N=taille_x*taille_y;
  fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
  fftw_plan p;
  //Réservation memoire
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  //Récupération de l'image dans la partie reelle de l'entree
  for(int cpt=0;cpt<N;cpt++)
    {
      in[cpt][0]=entree_reelle[cpt];
      in[cpt][1]=entree_imag[cpt];
    }
  //calcul du plan, parametre servant a calculer et optimiser le FFT
  p=fftw_plan_dft_2d( taille_x,  taille_y, in, out,FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p); /* repeat as needed */

  for(int cpt=0;cpt<(N);cpt++)
    {

      fft_reel[cpt]=out[cpt][0]/(N*N); //division par N^2 pour normaliser la fftw qui n'est pas normalisée
      fft_imag[cpt]=out[cpt][1]/(N*N);
    }
  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}


void TF2D_1P(double entree[],double sortie[],int taille_x,int taille_y)
{
  fftw_plan_with_nthreads(6);

  int N=taille_x*taille_y;
  fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
  fftw_plan p;
  //Réservation memoire
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  //Récupération de l'image dans la partie reelle de l'entree
  for(int cpt=0;cpt<N;cpt++)
    {
      in[cpt][0]=entree[cpt];
      in[cpt][1]=0.0f;
    }
  //calcul du plan, parametre servant a calculer et optimiser le FFT
  p=fftw_plan_dft_2d( taille_x,  taille_y, in, out,FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p); /* repeat as needed */

  for(int cpt=0;cpt<(N);cpt++)
    {
      sortie[cpt]=out[cpt][0]/N; //division par N^2 pour normaliser la fftw qui n'est pas normalisée
    }
  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}

void TF2D_INV(double entree_reelle[],double entree_imag[],double sortie_reelle[],double sortie_imag[],int taille_x,int taille_y)
{
  fftw_plan_with_nthreads(6);
  int N=taille_x*taille_y;
  fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
  fftw_plan p;
  //Réservation memoire
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  //Récupération de l'image dans la partie reelle de l'entree
  for(int cpt=0;cpt<N;cpt++)
    {
      in[cpt][0]=entree_reelle[cpt];
      in[cpt][1]=entree_imag[cpt];
    }
  //calcul du plan, parametre servant a calculer et optimiser le FFT
  p=fftw_plan_dft_2d(taille_x,  taille_y, in, out,FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p); /* repeat as needed */

  for(int cpt=0;cpt<(N);cpt++)
    {

      sortie_reelle[cpt]=out[cpt][0];
      sortie_imag[cpt]=out[cpt][1];
    }
  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}

void TF2D_INV_1P(double entree[],double sortie[],int taille_x,int taille_y)
{
  fftw_plan_with_nthreads(6);
  int N=taille_x*taille_y;
  fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
  fftw_plan p;
  //Réservation memoire
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  //Récupération de l'image dans la partie reelle de l'entree
  for(int cpt=0;cpt<N;cpt++)
    {
      in[cpt][0]=entree[cpt];
      in[cpt][1]=0.0f;
    }
  //calcul du plan, parametre servant a calculer et optimiser le FFT
  p=fftw_plan_dft_2d( taille_x,  taille_y, in, out,FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p); /* repeat as needed */

  for(int cpt=0;cpt<(N);cpt++)
    {
      sortie[cpt]=out[cpt][0];
    }
  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}


/////////////////////////////////////////////////////////////////////////////////////
//fonction de circshift

/////////////////////////////////////////////////////////////////////////////////////
void circshift(double* entree,double* entree_shift, int dimx,int decal_x)
{
  //int dimy=dimx;
  int decal_y=decal_x;

  for(int yi=0;yi<decal_y;yi++)
    {	for(int xi=0;xi<decal_x;xi++)

	{
	  int pixel=yi*dimx+xi;
	  int pixel_shift=(yi+decal_y)*dimx+xi+decal_x;
	  //1er quadrant vers 4 eme
	  entree_shift[pixel_shift]=entree[pixel];
	  //4 eme quadrant vers 1er
	  entree_shift[pixel]=entree[pixel_shift];
	  //2eme vers 3eme
	  entree_shift[(yi+decal_y)*dimx+xi]=entree[pixel+decal_x];
	  //3eme vers 2eme
	  entree_shift[pixel+decal_x]=entree[(yi+decal_y)*dimx+xi];
	}
    }
}

void circshift2(double* entree, double* result, Var2D dim,Var2D decal)
{

  //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
  decal.y=decal.y%dim.y;
  decal.x=decal.x%dim.x;

  for(int yi=0;yi<dim.y;yi++)
    {

      int yi_decal=yi+decal.y;
      //si dépassement des dimensions
      if(yi+decal.y>dim.y-1)
	yi_decal=-dim.y+yi+decal.y;
      //si décage négatif
      if(yi+decal.y<0)
	yi_decal=dim.y+yi+decal.y;
      int hauteur_decal=yi_decal*dim.x;
      int hauteur=yi*dim.x;
      for(int xi=0;xi<dim.x;xi++)
	{
	  int xi_decal=xi+decal.x;

	  if(xi+decal.x>dim.x-1)
	    xi_decal=-dim.x+xi+decal.x;
	  if(xi+decal.x<0)
	    xi_decal=dim.x+xi+decal.x;

	  int pos1D_init=hauteur+xi;
	  int pos1D_shift=hauteur_decal+xi_decal;//nouvelle position après décalage
	  result[pos1D_shift]=entree[pos1D_init];

	}
    }
}

void circshift3D2(double *volume3D, double *volume3D_shift, Var3D dimFinal3D, Var3D decal3D)
{
  decal3D.x=decal3D.x%dimFinal3D.x;//élmiiner les "modulos"
  decal3D.y=decal3D.y%dimFinal3D.y;
  decal3D.z=decal3D.z%dimFinal3D.y;

  unsigned short int xi,yi,zi=0;
  short int x2,y2,z2=0; //signé car une fois décalé, peuvent être négatifs!
  const unsigned int taille_plan =dimFinal3D.x*dimFinal3D.y;

  for(zi=0;zi<dimFinal3D.z;zi++)
    {
      if(zi+decal3D.z>dimFinal3D.z-1)//dépassement à droite
	{
	  z2=zi+decal3D.z-dimFinal3D.z;
	}
      else
	{
	  if(zi+decal3D.z<0)//dépassement à gauche
	    {
	      z2=dimFinal3D.z+(decal3D.z+zi);
	    }
	  else
	    {
	      z2=zi+decal3D.z;
	    }
	}
      int nb_pixelz_decal=z2*taille_plan;
      unsigned int nb_pixelz=zi*taille_plan;
      for(yi=0;yi<dimFinal3D.y;yi++)
	{
	  if(yi+decal3D.y>dimFinal3D.y-1)//dépassement à droite
	    {
	      y2=yi+decal3D.y-dimFinal3D.y;
	    }
	  else
	    {
	      if(yi+decal3D.y<0)//dépassement à gauche
		{
		  y2=dimFinal3D.y+(decal3D.y+yi);
		}
	      else
		{
		  y2=yi+decal3D.y;
		}
	    }
	  int nb_lignes=yi*dimFinal3D.x;
	  int nb_lignes_decal=y2*dimFinal3D.x;

	  for(xi=0;xi<dimFinal3D.x;xi++)
	    {
	      if(xi+decal3D.x>dimFinal3D.x-1)//dépassement à droite
                {
                  x2=xi+decal3D.x-dimFinal3D.x;
		}
	      else
		{
		  if(xi+decal3D.x<0)//dépassement à gauche
		    {
		      x2=dimFinal3D.x+(decal3D.x+xi);
		    }
		  else
		    {
		      x2=xi+decal3D.x;
		    }
		}

	      // pixel=zi*dimFinal3D.x*dimFinal3D.y+yi*dimFinal3D.x+xi;
	      //pixel_decal=z2*dimFinal3D.x*dimFinal3D.y+y2*dimFinal3D.x+x2;
	      //volume3D_shift[pixel_decal]=volume3D[pixel];
	      //volume3D_shift[z2*dimFinal3D.x*dimFinal3D.y+y2*dimFinal3D.x+x2]=volume3D[zi*dimFinal3D.x*dimFinal3D.y+yi*dimFinal3D.x+xi];
	      volume3D_shift[nb_pixelz_decal+nb_lignes_decal+x2]=volume3D[nb_pixelz+nb_lignes+xi];
	      //memcpy(volume3D_shift+(nb_pixelz_decal+nb_lignes_decal+x2), volume3D+(nb_pixelz+nb_lignes+xi), 8); commande ok mais plus lent !?

	    }
	}
    }
}
///fonction de masquage
void multiplier_masque(double image[], unsigned char masque[], int t_image, int t_mask, int centreX, int centreY)
{int t_imageX=t_image;
  int t_imageY=t_image;
  //if((t_mask%2)!=0)
  //cout<<"fonction multiplier_masque : attention, la dimension de votre masque est impaire! t_mask="<<t_mask<<endl;

  for(int x=centreX-t_mask/2;x<centreX+t_mask/2;x++)
    {
      for(int y=centreY-t_mask/2;y<centreY+t_mask/2;y++)
	{   //si le masque déborde de l'image (attention le masque repasse du coté gauche lorsque le jumeau touche le bord droit! A corriger
	  int x_2=x;
	  int y_2=y;
	  if(x>2*t_imageX)
	    x_2=x-2*t_imageX;
	  if(y>2*t_imageY)
	    y_2=y-2*t_imageY;
	  if(x<0)
	    x_2=2*t_imageX+x;
	  if(y<0)
	    y_2=2*t_imageY+y;
	  ///coordonnées 1D du centre
	  int cptj=2*t_imageY*y_2+x_2;
	  //if((double(masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)])/255)!=1)
	  ///cout<<"masque:"<<(double(masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)])/255)<<endl;
	  image[cptj]=image[cptj]*masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)]/255;
	}
    }
}

///fonction de masquage
void multiplier_masque2(double image[], double masque[], int t_image, int t_mask, int centreX, int centreY)
{int t_imageX=t_image;
  int t_imageY=t_image;
  //if((t_mask%2)!=0)
  //cout<<"fonction multiplier_masque : attention, la dimension de votre masque est impaire! t_mask="<<t_mask<<endl;

  for(int x=centreX-t_mask/2;x<centreX+t_mask/2;x++)
    {
      for(int y=centreY-t_mask/2;y<centreY+t_mask/2;y++)
	{   //si le masque déborde de l'image (attention le masque repasse du coté gauche lorsque le jumeau touche le bord droit! A corriger
	  int x_2=x;
	  int y_2=y;
	  if(x>2*t_imageX)
	    x_2=x-2*t_imageX;
	  if(y>2*t_imageY)
	    y_2=y-2*t_imageY;
	  if(x<0)
	    x_2=2*t_imageX+x;
	  if(y<0)
	    y_2=2*t_imageY+y;
	  ///coordonnées 1D du centre
	  int cptj=2*t_imageY*y_2+x_2;
	  //if((double(masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)])/255)!=1)
	  ///cout<<"masque:"<<(double(masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)])/255)<<endl;
	  image[cptj]=image[cptj]*masque[t_mask*(y-centreY+t_mask/2)+(x-centreX+t_mask/2)];
	}
    }
}
void crop_window(double* img_src, size_t src_dimx, size_t src_dimy, double* img_dst, size_t dst_dimx, size_t dst_dimy, size_t ul_corner_x, size_t ul_corner_y)
{
  assert(dst_dimx < src_dimx);
  assert(dst_dimy < src_dimy);

  assert((ul_corner_x + dst_dimx) <= src_dimx);
  assert((ul_corner_y + dst_dimy) <= src_dimy);

  double *start_src = img_src + (ul_corner_y * src_dimx);

  double *curs_src = start_src;
  size_t pos_dst = 0;
  size_t i, j;

  for (i = 0; i < dst_dimy; i++) {
    for (j = ul_corner_x; j < ul_corner_x + dst_dimx; j++) {
      img_dst[pos_dst] = curs_src[j];
      pos_dst++;
    }
    curs_src += src_dimx;
  }
}
///phase2pi
void phase2pi(double* phaseMod1pi, double* partie_reel, double* partie_imag, int taille_x, int taille_y)
{

  ///calcul phase de -PI à PI avec argTan, attention 0 en phase
  /*
    for(int pixel=0;pixel<taille_x*taille_y;pixel++)
    {
    double denominPhase=partie_reel[pixel];
    double numeraPhase=partie_imag[pixel];

    if(denominPhase!=0)
    {
    double argTanPhase=numeraPhase/denominPhase;
    phaseMod1pi[pixel]=atan(argTanPhase);

    }
    else
    {
    phaseMod1pi[pixel]=0;
    }
    }
  */
  ///calcul phase de -PI à PI avec argSin
  for(int pixel=0;pixel<taille_x*taille_y;pixel++)
    {
      double denominPhase=partie_reel[pixel];
      double numeraPhase=partie_imag[pixel];

      if(denominPhase>0)
	{
	  double racine=sqrt(pow(numeraPhase,2)+pow(denominPhase,2));
	  phaseMod1pi[pixel]=asin(numeraPhase/racine);
	}
      else
	{
	  double racine=sqrt(pow(numeraPhase,2)+pow(denominPhase,2));
	  phaseMod1pi[pixel]=PI-asin(numeraPhase/racine);
	}
    }

  ///calcul phase de -pi à PI
  //     for(int pixel=0;pixel<taille_x*taille_y;pixel++)
  //     {
  //        double cos_phase=partie_reel[pixel];
  //        double sin_phase=partie_imag[pixel];
  //        if(sin_phase>0)
  //        {
  //            if(cos_phase>0)
  //            {
  //                double argTanPhase=sin_phase/cos_phase;
  //                phaseMod1pi[pixel]=atan(argTanPhase);
  //            }
  //            else if(cos_phase<0)
  //            {
  //                double argTanPhase=sin_phase/cos_phase;
  //                phaseMod1pi[pixel]=atan(argTanPhase)+PI;
  //            }
  //            else if(cos_phase==0)
  //            {
  //                phaseMod1pi[pixel]=PI/2;
  //            }
  //        }
  //        else if(sin_phase<0)
  //        {
  //            if(cos_phase<0)
  //            {
  //                double argTanPhase=sin_phase/cos_phase;
  //                phaseMod1pi[pixel]=atan(argTanPhase)-PI;
  //            }
  //            else if(cos_phase>0)
  //            {
  //                double argTanPhase=sin_phase/cos_phase;
  //                phaseMod1pi[pixel]=atan(argTanPhase);
  //            }
  //            else if(cos_phase==0)
  //            {
  //                phaseMod1pi[pixel]=3*PI/2;
  //            }
  //        }
  //        else if(sin_phase==0)
  //        {
  //            if(cos_phase>0)
  //            {
  //                phaseMod1pi[pixel]=0;
  //            }
  //            else if(cos_phase<0)
  //            {
  //                phaseMod1pi[pixel]=PI;
  //            }
  //
  //        }
  //     }

}

///phase2pi
void phaseUnwrapping(double* partie_reel, double* partie_imag, double* WrappedImage,
                     double* UnwrappedImage, PIXEL *Pixel, EDGE *Edge, int taille)
{
  int image_height=taille, image_width=taille;
  int image_size = image_height * image_width;
  //int two_image_size = 2 * image_size;
  int No_of_Edges = (image_width)*(image_height-1) + (image_width-1)*(image_height);

  phase2pi(WrappedImage, partie_reel, partie_imag, taille, taille);

  //initialise the pixels
  initialisePIXELs(WrappedImage, Pixel, image_width, image_height);

  calculate_reliability(WrappedImage, Pixel, image_width, image_height);

  horizentalEDGEs(Pixel, Edge, image_width, image_height);

  verticalEDGEs(Pixel, Edge, image_width, image_height);

  //sort the EDGEs depending on their reiability. The PIXELs with higher relibility (small value) first
  //if your code stuck because of the quicker_sort() function, then use the quick_sort() function
  //run only one of the two functions (quick_sort() or quicker_sort() )
  quick_sort(Edge, No_of_Edges);
  //quicker_sort(edge, edge + No_of_Edges - 1);

  //gather PIXELs into groups
  gatherPIXELs(Edge, image_width, image_height);

  //unwrap the whole image
  unwrapImage(Pixel, image_width, image_height);
  //copy the image from PIXEL structure to the wrapped phase array passed to this function
  returnImage(Pixel, UnwrappedImage, image_width, image_height);
}

///#### la correction de fond après élimination des franges d'inclinaision avec ajustement polynomiale ######
void corrplan_I(double *plan_reel, double *plan_imag, double *plan_reel_fit, double *plan_imag_fit, int taille,  Mat masqueC)
{
  ///Copier double vers matrice de opencv///
  Mat imgR(taille,taille,CV_64F);
  Mat imgI(taille,taille,CV_64F);
  for(int y=0; y<taille; y++)
    {
      for(int x=0; x<taille; x++)
	{
	  int cpt=y*taille+x;
	  imgR.at<double>(y,x)=plan_reel[cpt];
	  imgI.at<double>(y,x)=plan_imag[cpt];
	}
    }

  int degpoly=3;
  Mat CoefsolveR;
  Mat CoefsolveI;

  CalculPoly(imgR,masqueC,CoefsolveR,degpoly,true);
  CalculPoly(imgI,masqueC,CoefsolveI,degpoly,true);

  Mat resultatpolyR(imgR.rows,imgR.cols,CV_64F);
  Mat resultatpolyI(imgI.rows,imgI.cols,CV_64F);

  CalculFond(CoefsolveR,degpoly,resultatpolyR);
  CalculFond(CoefsolveI,degpoly,resultatpolyI);
  write_dataM("/ramdisk/fond_reel.bin", resultatpolyR);

  Mat tmpR, tmpI;
  imgR.convertTo(tmpR, CV_64F);
  imgI.convertTo(tmpI, CV_64F);

  Mat plan_reelFit(imgR.rows,imgR.cols,CV_64F);
  Mat plan_imagFit(imgI.rows,imgI.cols,CV_64F);

  divComplx(tmpR, tmpI, resultatpolyR, resultatpolyI, plan_reelFit, plan_imagFit);
  //write_dataM("/home/hui/maniptomo/IDP/champ_2d/plan_reel_fit.bin", plan_reelFit);

  ///Copier matrice de opencv vers double///
  for(int y=0; y<plan_reelFit.rows; y++)
    {
      for(int x=0; x<plan_reelFit.cols; x++)
	{
	  int cpt=y*plan_reelFit.cols+x;
	  plan_reel_fit[cpt]=plan_reelFit.at<double>(y,x);
	  plan_imag_fit[cpt]=plan_imagFit.at<double>(y,x);
	}
    }
}

void corrPhase(double *phase, double *phase_fit, int taille,  Mat masqueC)
{
  ///Copier double vers matrice de opencv///
  Mat img(taille,taille,CV_64F);
  for(int y=0; y<taille; y++)
    {
      for(int x=0; x<taille; x++)
	{
	  int cpt=y*taille+x;
	  img.at<double>(y,x)=phase[cpt];

	}
    }

  int degpoly=3;
  Mat Coefsolve;
  CalculPoly(img,masqueC,Coefsolve,degpoly,true);
  //cout<<Coefsolve<<endl;
  Mat resultatpoly(img.rows,img.cols,CV_64F);
  CalculFond(Coefsolve,degpoly,resultatpoly);
  write_dataM("/ramdisk/fond_reel.bin", resultatpoly);

  Mat tmp;
  img.convertTo(tmp, CV_64F);
  tmp=tmp-resultatpoly;

  ///Copier matrice de opencv vers double///
  for(int y=0; y<tmp.rows; y++)
    {
      for(int x=0; x<tmp.cols; x++)
	{
	  int cpt=y*tmp.cols+x;
	  phase_fit[cpt]=tmp.at<double>(y,x);
	}
    }
}

void corrAmpli(double *ampli, double *ampli_fit, int taille,  Mat masqueC)
{
  ///Copier double vers matrice de opencv///
  Mat img(taille,taille,CV_64F);
  for(int y=0; y<taille; y++)
    {
      for(int x=0; x<taille; x++)
	{
	  int cpt=y*taille+x;
	  img.at<double>(y,x)=ampli[cpt];

	}
    }

  int degpoly=3;
  Mat Coefsolve;
  CalculPoly(img,masqueC,Coefsolve,degpoly,true);
  Mat resultatpoly(img.rows,img.cols,CV_64F);
  CalculFond(Coefsolve,degpoly,resultatpoly);
  //write_dataM("/home/hui/maniptomo/IDP/champ_2d/fond_reel.bin", resultatpoly);

  Mat tmp;
  img.convertTo(tmp, CV_64F);
  tmp=tmp/resultatpoly;

  ///Copier matrice de opencv vers double///
  for(int y=0; y<tmp.rows; y++)
    {
      for(int x=0; x<tmp.cols; x++)
	{
	  int cpt=y*tmp.cols+x;
	  ampli_fit[cpt]=tmp.at<double>(y,x);
	}
    }
}

///########################################fonctions en dessous sertent à corriger le fond#################################################
double Poly2DEval(Mat Coefficients, int deg,int x, int y)
{
  int j=0,k=0;
  double somme=0;
  for(int i=0; i<=deg; i++)
    {
      j=0;

      while((i+j)<=deg)
        {
	  somme+=(Coefficients.at<double>(k))*pow(x,i)*pow(y,j);
	  k++;
	  j++;
        }

    }
  return somme;
}

int TaillePoly2D(int deg)
{
  int j=0;
  int taille=0;
  for(int i=0; i<=deg; i++)
    {

      while((i+j)<=deg)
        {
	  taille++;
	  j++;
        }
      j=0;
    }
  return taille;
}

int CountM(Mat& Masque)
{
  int count=0;
  for(int i=0; i<Masque.rows; i++)
    {
      for(int j=0; j<Masque.cols; j++)
        {
	  if(Masque.at<uchar>(j,i)>50)
            count++;
        }
    }
  //cout<<"count="<<count<<endl;
  return count;
}

void CalculPoly(Mat& imagebrute,Mat& Masque, Mat& polynome,int deg, bool method)
{

  int nbcolonnes=TaillePoly2D(deg);
  //int nblignes= countNonZero(Masque);
  int nblignes= CountM(Masque);
  //int row=0;
  int col=0;
  int k = 0;
  Mat B(Size(nbcolonnes,nblignes),CV_64F);
  Mat Bt(Size(nblignes,nbcolonnes),CV_64F);
  Mat f(nblignes,1,CV_64F);

  for (int y=0; y< imagebrute.rows; y++)
    {
      for (int x=0; x<imagebrute.cols; x++)
        {
	  col=0;
	  if (Masque.at<uchar>(y,x)>50)
            {

	      f.at<double>(k)=(double)imagebrute.at<double>(y,x);

	      for(int i=0; i<=deg; i++)
                {
		  int j=0;
		  while((i+j)<=deg)
                    {
		      B.at<double>(k,col)=pow((double)x,i)*pow((double)y,j);
		      col++;

		      j++;
                    }
                }
	      k++;
            }
        }
    }

  Mat coef(nbcolonnes,1,CV_64F);
  Mat D(Size(nbcolonnes,nblignes),CV_64F);
  Mat invD(Size(nbcolonnes,nblignes),CV_64F);

  if (method)  cv::solve(B,f,coef,DECOMP_NORMAL);
  else
    {
      cv::transpose(B,Bt);
      D=Bt*B;
      cv::invert(D,invD);
      coef= (invD*Bt)*f;
    }
  coef.copyTo(polynome);
}

void CalculFond(Mat coefficients, int deg, Mat imagefond)

{
  int col=0;
  int j=0;
  double somme=0;

  for (int y=0; y<imagefond.rows; y++)
    {
      for (int x=0; x<imagefond.cols; x++)
        {
	  col=0;
	  somme=0;

	  for(int i=0; i<=deg; i++)
            {
	      j=0;
	      while((i+j)<=deg)
                {
		  somme+= (coefficients.at<double>(col))*pow((double)x,i)*pow((double)y,j);
		  col++;
		  j++;
                }
            }

	  imagefond.at<double>(y,x)=Poly2DEval(coefficients,deg,x,y);

        }
    }
}

void read_dataM(char *inputfile,Mat* Data, int image_size, int IMAGE_NUM)
{
  //printf("Reading the Data Values form Binary File.............>");
  FILE *ifptr;
  float precision;
  Mat Data_tmp(image_size,image_size, CV_64F);
  ifptr = fopen(inputfile,"rb");
  if(ifptr == NULL) printf("Error opening the file\n");
  for(int z=0;z<IMAGE_NUM;z++)
    {
      for(int i=0;i<image_size;i++)
        {
	  for(int j=0; j< image_size; j++)
            {
	      fread(&precision,sizeof(precision),1,ifptr);
	      Data_tmp.at<double>(i,j)=precision;
            }
        }
      //cout<<"precision:"<<z<<endl;
      Data_tmp.copyTo(Data[z]);
    }
  fclose(ifptr);
  //printf(" Done.\n");
}

void write_dataM(char *outputfile,Mat Data)
{
  //printf("Writing the Unwrapped Values to Binary File.............>");
  FILE *ifptr;
  ifptr = fopen(outputfile,"a+b");
  if(ifptr == NULL) printf("Error opening the file\n");

  for(int i=0;i<Data.rows;i++)
    for(int j=0; j< Data.cols; j++)
      {
        {
	  float precision=(float)Data.at<double>(i,j);
	  fwrite(&precision,sizeof(precision),1,ifptr);
        }
      }

  fclose(ifptr);
  //printf(" Done.\n");
}

void divComplx(Mat a, Mat b, Mat c, Mat d, Mat& x, Mat& y)
{
  Mat mod, c_pow, d_pow;
  pow(c,2,c_pow);
  pow(d,2,d_pow);
  mod = c_pow + d_pow;
  x= (a.mul(c) + b.mul(d))/mod;
  y= (b.mul(c) - a.mul(d))/mod;
}

Var2D corr_crois(double* obj2D_A, double* obj2D_B, Var2D dim)
{
  int NbPts2D=dim.x*dim.y;

  //copie pour passer en complexe...
  nbCplx *copie_obj2D_A=new nbCplx[NbPts2D];
  nbCplx *copie_obj2D_B=new nbCplx[NbPts2D];
  for(int cpt=0;cpt<NbPts2D;cpt++)//copies pour passer en complexe (la fonction Tf2D demandent des cplx)
    {
      copie_obj2D_A[cpt].Re=obj2D_A[cpt];
      copie_obj2D_A[cpt].Im=0;
      copie_obj2D_B[cpt].Re=obj2D_B[cpt];
      copie_obj2D_B[cpt].Im=0;
    }
  //stockage des spectres
  nbCplx *spect2D_A= new nbCplx[NbPts2D];
  nbCplx *spect2D_B= new nbCplx[NbPts2D];
  //stocker la correlation et son spectre
  nbCplx *spectCorr=new nbCplx[NbPts2D];
  nbCplx *Corr=new nbCplx[NbPts2D];


  //SAV_Re(copie_obj2D_B,NbPts2D,"/home/mat/tomo_test/copie_objetB.bin",FLOAT,"wb");
  TF2Dcplx(copie_obj2D_A,spect2D_A,dim);
  TF2Dcplx(copie_obj2D_B,spect2D_B,dim);

  //Elimination du module pour isoler  la phase : A*conj(B)/(modA*modB)

  for(int cpt=0;cpt<NbPts2D;cpt++){
    //double modA=sqrt(pow(spect2D_A[cpt].Re,2)+pow(spect2D_A[cpt].Im,2));
    //double modB=sqrt(pow(spect2D_B[cpt].Re,2)+pow(spect2D_B[cpt].Im,2));
    spectCorr[cpt].Re=double(spect2D_A[cpt].Re*spect2D_B[cpt].Re+spect2D_A[cpt].Im*spect2D_B[cpt].Im);
    spectCorr[cpt].Im=double(spect2D_A[cpt].Im*spect2D_B[cpt].Re-spect2D_A[cpt].Re*spect2D_B[cpt].Im);
  }
  //TF de l'expoentielle contenant le dépahsage->translation dans l'espace direct
  TF2Dcplx_INV(spectCorr,Corr,dim);
  int cptMax=0;
  double valMax=0;
  double valModCorr=0;
  for(int cpt=0;cpt<NbPts2D;cpt++)
    {
      valModCorr=pow(Corr[cpt].Re,2)+pow(Corr[cpt].Im,2);

      if(valModCorr>valMax)
        {
	  valMax=valModCorr;
	  cptMax=cpt;
        }
    }
  Var2D decal2D_I={0,0};
  decal2D_I.x=cptMax%dim.x;
  decal2D_I.y=cptMax/dim.y;
  if(decal2D_I.x>15 && decal2D_I.x<dim.x-15)
    decal2D_I.x=0;
  if(decal2D_I.y>15 &&  decal2D_I.y<dim.y-15)
    decal2D_I.y=0;
  //cout<<"---------------------------"<<endl;
  //cout<<"X="<<decal2D_I.x<<endl;
  //cout<<"Y="<<decal2D_I.y<<endl;
  // SAV_Re(spectCorr,NbPts2D,"/home/mat/tomo_test/SpectCorrRE.bin",FLOAT,"wb");
  //SAV_Re(Corr,NbPts2D,"/home/mat/tomo_test/Corr.bin",FLOAT,"wb");
  return decal2D_I;
}

void decal2DGen(double* entree, double* result, Var2D dim,Var2D decal)
{
  //si décalage supérieure à dim, on fait plus d'un tour, donc on prend le modulo
  decal.y=decal.y%dim.y;
  decal.x=decal.x%dim.x;
  if(decal.x<0)
    decal.x=dim.x+decal.x;
  if(decal.y<0)
    decal.y=dim.y+decal.y;
  for(int yi=0; yi<dim.y-decal.y; yi++) {
    for(int xi=0; xi<dim.x-decal.x; xi++)
      { //cout<<"xi,yi="<<xi<<","<<yi<<endl;
	int pixel=yi*dim.x+xi;
	int pixel_shift=(yi+decal.y)*dim.x+xi+decal.x;
	result[pixel_shift]=entree[pixel];
      }
    for(int xi=dim.x-decal.x; xi<dim.x; xi++)
      {
	int pixel=yi*dim.x+xi;
	int pixel_shift=(yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
	result[pixel_shift]=entree[pixel];
      }
  }
  for(int yi=dim.y-decal.y; yi<dim.y; yi++) {
    for(int xi=0; xi<dim.x-decal.x; xi++)
      {
	int pixel=yi*dim.x+xi;
	int pixel_shift=(-dim.y+yi+decal.y)*dim.x+xi+decal.x;
	result[pixel_shift]=entree[pixel];
      }
    for(int xi=dim.x-decal.x; xi<dim.x; xi++)
      {
	int pixel=yi*dim.x+xi;
	int pixel_shift=(-dim.y+yi+decal.y)*dim.x+(-dim.x+xi+decal.x);//Y_shift*dim.x+X_shift;
	result[pixel_shift]=entree[pixel];
      }
  }

}
void TF2Dcplx(nbCplx *entree, nbCplx *fft, Var2D dim)
{
  fftw_plan_with_nthreads(3);
  int N=dim.x*dim.y;
  fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
  fftw_plan p;
  //Réservation memoire
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  //Récupération de l'image dans la partie reelle de l'entree
  for(int cpt=0; cpt<N; cpt++) {
    in[cpt][0]=entree[cpt].Re;
    in[cpt][1]=entree[cpt].Im;
  }
  //calcul du plan, parametre servant a calculer et optimiser le FFT
  p=fftw_plan_dft_2d( dim.x,  dim.y, in, out,FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p); /* repeat as needed */

  for(int cpt=0; cpt<(N); cpt++) {
    fft[cpt].Re=out[cpt][0]; //division par N^2 pour normaliser la fftw qui n'est pas normalisée
    fft[cpt].Im=out[cpt][1];
  }
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
}
void TF2Dcplx_INV(nbCplx *fft_entree, nbCplx *objet, Var2D dim)
{
  fftw_plan_with_nthreads(3);
  int N=dim.x*dim.y;
  fftw_complex *in, *out;//Déclaration des variables pour la FFT : entree,sortie et "fftplan"
  fftw_plan p;
  //Réservation memoire
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  //Récupération de l'image dans la partie reelle de l'entree
  for(int cpt=0; cpt<N; cpt++) {
    in[cpt][0]=fft_entree[cpt].Re;
    in[cpt][1]=fft_entree[cpt].Im;
  }
  //calcul du plan, parametre servant a calculer et optimiser le FFT
  p=fftw_plan_dft_2d( dim.x,  dim.y, in, out,FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p); /* repeat as needed */

  for(int cpt=0; cpt<(N); cpt++) {
    objet[cpt].Re=out[cpt][0]; //division par N^2 pour normaliser la fftw qui n'est pas normalisée
    objet[cpt].Im=out[cpt][1];
  }
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
}

void Resharp(double* entree, double* sortie, int taille)
{
  Mat imgResult;
  ///Copier double vers matrice de opencv///
  Mat img(taille,taille,CV_64F);
  for(int y=0; y<taille; y++)
    {
      for(int x=0; x<taille; x++)
	{
	  int cpt=y*taille+x;
	  img.at<double>(y,x)=entree[cpt];

	}
    }


  //    /// laplacian
  //    Mat imgLaplacian, imgResult;
  //
  //    // ok, now try different kernel
  //    Mat kernel = (Mat_<double>(3,3) <<
  ////        1,  1, 1,
  ////        1, -8, 1,
  ////        1,  1, 1);
  //        0,  1, 0,
  //        1, -4, 1,
  //        0,  1, 0);
  //
  //    filter2D(img, imgLaplacian, CV_64F, kernel);
  //    imgResult = img - imgLaplacian;


  //    /// sobel
  //    int scale = 1;
  //    int delta = 0;
  //    int ddepth = CV_64F;
  //
  //    Mat grad_x, grad_y;
  //
  //    /// Gradient X
  //    //Scharr( img1[0], grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
  //    Sobel( img, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
  //    /// Gradient Y
  //    //Scharr( img1[0], grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
  //    Sobel( img, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
  //    /// Total Gradient (approximate)
  //    addWeighted( grad_x, 0.5, grad_y, 0.5, 0, imgResult );

  /// laplacian
  int kernel_size = 3;
  int scale = 1;
  int delta = 0;
  int ddepth = CV_64F;
  Laplacian(img, imgResult, ddepth, kernel_size, scale, delta, BORDER_DEFAULT );

  ///Copier matrice de opencv vers double///
  for(int y=0; y<imgResult.rows; y++)
    {
      for(int x=0; x<imgResult.cols; x++)
	{
	  int cpt=y*imgResult.cols+x;
	  sortie[cpt]=imgResult.at<double>(y,x);
	}
    }
}

void recalAmpPha(double* amplitude, double* phase, double * ampliRef, double * phaseRef, Var2D dimFin, double* sortie_reel, double* sortie_imag,
                 int cpt_angle, Var2D posDecalAmpli, Var2D posDecalPhase, int NbAngle)
{
  ///#######################recal phase#####################///

  double *phase_recal= new double[dimFin.x*dimFin.y];

  if(cpt_angle==0)
    {
      for(int cpt=0;cpt<dimFin.x*dimFin.y;cpt++)
        {
	  phaseRef[cpt]=phase[cpt];
        }
    }
  posDecalPhase=corr_crois(phaseRef, phase, dimFin);
  decal2DGen(phase, phase_recal,dimFin, posDecalPhase);
  SAV(phase_recal, dimFin.x*dimFin.y, "/ramdisk/phaseMod1pi_recal.bin", FLOAT,"a+b");

  ///#######################recal amplitude#####################///

  double *ampli_recal= new double[dimFin.x*dimFin.y];

  if(cpt_angle==0)
    {
      for(int cpt=0;cpt<dimFin.x*dimFin.y;cpt++)
        {
	  ampliRef[cpt]=amplitude[cpt];
        }
    }
  posDecalAmpli=corr_crois(ampliRef, amplitude, dimFin);
  decal2DGen(amplitude, ampli_recal,dimFin, posDecalAmpli);
  SAV(ampli_recal, dimFin.x*dimFin.y, "/ramdisk/amplitude2D_recal.bin", FLOAT,"a+b");

  ///#######################recréer plan reel et imaginaire à partir de amplitude et phase recalées#####################///

  double *plan_reel_recal= new double[dimFin.x*dimFin.y];
  double *plan_imag_recal= new double[dimFin.x*dimFin.y];

  for(int cpt=0;cpt<dimFin.x*dimFin.y;cpt++)
    {
      plan_reel_recal[cpt]=ampli_recal[cpt]*cos(phase_recal[cpt]);
      plan_imag_recal[cpt]=ampli_recal[cpt]*sin(phase_recal[cpt]);
    }
  SAV(plan_reel_recal, dimFin.x*dimFin.y, "/ramdisk/plan_reel_recal.bin", FLOAT,"a+b");

  double *spectre_reel= new double[dimFin.x*dimFin.y];
  double *spectre_imag= new double[dimFin.x*dimFin.y];

  TF2D(plan_reel_recal,plan_imag_recal,spectre_reel,spectre_imag,dimFin.x,dimFin.y);

  int cpt_maxRecal=0;
  double *spectre_module= new double[dimFin.x*dimFin.y];
  for(int cpt=0;cpt<(dimFin.x*dimFin.y);cpt++)
    {
      spectre_module[cpt]=pow(spectre_reel[cpt],2)+pow(spectre_imag[cpt],2);
      if(spectre_module[cpt]>spectre_module[cpt_maxRecal])
        {
	  cpt_maxRecal=cpt;
        }
    }

  double max_reel = spectre_reel[cpt_maxRecal];
  double max_imag = spectre_imag[cpt_maxRecal];
  double max_modul = spectre_module[cpt_maxRecal];

  double *spectre_reel_norm= new double[dimFin.x*dimFin.y];
  double *spectre_imag_norm= new double[dimFin.x*dimFin.y];

  for(int cpt=0;cpt<(dimFin.x*dimFin.y);cpt++)
    {
      spectre_reel_norm[cpt]=(spectre_reel[cpt]*max_reel+spectre_imag[cpt]*max_imag)/max_modul;
      spectre_imag_norm[cpt]=(spectre_imag[cpt]*max_reel-spectre_reel[cpt]*max_imag)/max_modul;
      //printf("fft_reel_shift_norm: %f, fft_imag_shift_norm: %f \n ",fft_reel_shift_norm[cpt], fft_imag_shift_norm[cpt]);
    }

  TF2D_INV(spectre_reel_norm, spectre_imag_norm, sortie_reel, sortie_imag, dimFin.x, dimFin.y);

  delete[] ampli_recal, phase_recal, plan_reel_recal, plan_imag_recal, spectre_reel, spectre_imag, spectre_module, spectre_reel_norm, spectre_imag_norm;
}

void calculPlan(unsigned char* holo1, unsigned char* holo2, unsigned char* holo3, unsigned char* holo4, double* plan_reel, double* plan_imag,
                int taille, Var2D coin, char* Fichier_holo, double* masque, unsigned char* mod_ref)
{
  holo1=charger_image2D(holo1,0,Fichier_holo, coin.x, coin.y,taille,taille);
  holo2=charger_image2D(holo2,1,Fichier_holo, coin.x, coin.y,taille,taille);
  holo3=charger_image2D(holo3,2,Fichier_holo, coin.x, coin.y,taille,taille);
  holo4=charger_image2D(holo4,3,Fichier_holo, coin.x, coin.y,taille,taille);

  double* ampli_ref=new double[taille*taille];
  for(int pixel=0;pixel<taille*taille;pixel++)
    {
      plan_reel[pixel]=((double)holo1[pixel]-(double)holo3[pixel])*(masque[pixel]);
      plan_imag[pixel]=((double)holo4[pixel]-(double)holo2[pixel])*(masque[pixel]);
      //ampli_ref[pixel]= sqrt((double)mod_ref[pixel]);
      //plan_reel[pixel]=(((double)holo1[pixel]-(double)holo3[pixel]))*(masque[pixel])/ampli_ref[pixel];
      //plan_imag[pixel]=(((double)holo4[pixel]-(double)holo2[pixel]))*(masque[pixel])/ampli_ref[pixel];
    }
  SAV(plan_reel, taille*taille, "/ramdisk/plan_reel.bin", FLOAT,"a+b");
  delete[] ampli_ref;
}

void calculPlan1(unsigned char* holo1, unsigned char* holo2, unsigned char* holo3, unsigned char* holo4, double* plan_reel, double* plan_imag,
		 int taille, Var2D coin, char* Fichier_holo, double* masque, unsigned char* mod_ref, double * supRedonPlanN, double * supRedonPlan, int dimfinal)
{
  int N=taille*taille;
  holo1=charger_image2D(holo1,0,Fichier_holo, coin.x, coin.y,taille,taille);
  holo2=charger_image2D(holo2,1,Fichier_holo, coin.x, coin.y,taille,taille);
  holo3=charger_image2D(holo3,2,Fichier_holo, coin.x, coin.y,taille,taille);
  holo4=charger_image2D(holo4,3,Fichier_holo, coin.x, coin.y,taille,taille);

  //ContrSpeckle(double(holo1), taille);

  double *DataModu=new double[N];
  double* ampli_ref=new double[N];

  for(int pixel=0;pixel<N;pixel++)
    {
      ///qualité des hologrammes
      DataModu[pixel]=2*sqrt((pow(holo4[pixel]-holo2[pixel],2)+pow(holo1[pixel]-holo3[pixel],2)))/(holo1[pixel]+holo2[pixel]+holo3[pixel]+holo4[pixel]);
      //cout<<DataModu[pixel]<<endl;
      if(DataModu[pixel]>0.1)
        {
	  plan_reel[pixel]=((double)holo1[pixel]-(double)holo3[pixel])*(masque[pixel]);
	  plan_imag[pixel]=((double)holo4[pixel]-(double)holo2[pixel])*(masque[pixel]);
	  //ampli_ref[pixel]= sqrt((double)mod_ref[pixel]);
	  //plan_reel[pixel]=(((double)holo1[pixel]-(double)holo3[pixel]))*(masque[pixel])/ampli_ref[pixel];
	  //plan_imag[pixel]=(((double)holo4[pixel]-(double)holo2[pixel]))*(masque[pixel])/ampli_ref[pixel];
	  supRedonPlanN[pixel]+=1;
        }
      else
        {
	  plan_reel[pixel]=0;
	  plan_imag[pixel]=0;

        }

    }
  SAV(DataModu, N, "/ramdisk/DataModu.bin", FLOAT,"a+b");
  SAV(plan_reel, N, "/ramdisk/plan_reel.bin", FLOAT,"a+b");

  //SAV(supRedonPlanN, N, "/home/hui/maniptomo/IDP/champ_2d/supRedonP.bin", FLOAT,"a+b");

  ///Copier double vers matrice de opencv///
  Mat img(taille,taille,CV_8U);
  for(int y=0; y<taille; y++)
    {
      for(int x=0; x<taille; x++)
	{
	  int cpt=y*taille+x;
	  img.at<char>(y,x)=supRedonPlanN[cpt];
	}
    }

  int size=dimfinal;
  Mat tmp(size,size,CV_8U);
  resize(img, tmp, tmp.size(), 0, 0);

  ///Copier matrice de opencv vers double///
  for(int y=0; y<tmp.rows; y++)
    {
      for(int x=0; x<tmp.cols; x++)
	{
	  int cpt=y*tmp.cols+x;
	  supRedonPlan[cpt]=tmp.at<char>(y,x);
	}
    }
  SAV(supRedonPlan, size*size, "/ramdisk/supRedonPlan.bin", FLOAT,"wb");
  delete[] ampli_ref;
}

void calculPlanCritLG(unsigned char* holo1, unsigned char* holo2, unsigned char* holo3, unsigned char* holo4, double* plan_reel, double* plan_imag,
		      int taille, Var2D coin, char* Fichier_holo, double* masque, unsigned char* mod_ref, int cpt_angle)
{
  int N=taille*taille;
  holo1=charger_image2D(holo1,0,Fichier_holo, coin.x, coin.y,taille,taille);
  holo2=charger_image2D(holo2,1,Fichier_holo, coin.x, coin.y,taille,taille);
  holo3=charger_image2D(holo3,2,Fichier_holo, coin.x, coin.y,taille,taille);
  holo4=charger_image2D(holo4,3,Fichier_holo, coin.x, coin.y,taille,taille);

  double *DataModu=new double[N];
  double *sommeP_reel=new double[N];
  double *sommeP_imag=new double[N];
  double *plan_I42=new double[N];
  double *plan_I13=new double[N];
  double *supRedonQ=new double[N];

  memset(plan_I42, 0, N*8);
  memset(plan_I13, 0, N*8);
  for(int pixel=0;pixel<N;pixel++)
    {
      plan_I13[pixel]=(double)holo1[pixel]-(double)holo3[pixel];
      plan_I42[pixel]=(double)holo4[pixel]-(double)holo2[pixel];
    }

  memset(DataModu, 0, N*8);
  memset(sommeP_reel, 0, N*8);
  memset(sommeP_imag, 0, N*8);

  //int pixelQ_max=0;
  double sommeQholo=0;
  double moyQholo=0;
  for(int pixel=0;pixel<N;pixel++)
    {
      ///qualité des hologrammes
      DataModu[pixel]=2*sqrt((pow(holo4[pixel]-holo2[pixel],2)+pow(holo1[pixel]-holo3[pixel],2)))/(holo1[pixel]+holo2[pixel]+holo3[pixel]+holo4[pixel]);
      //cout<<"contraste pour chaque pixel: "<<DataModu[pixel]<<endl;
      /*if(DataModu[pixel]>DataModu[pixelQ_max])
        {
	pixelQ_max=pixel;
        }*/
      if(DataModu[pixel]>0.0)
        {
	  sommeQholo=sommeQholo+DataModu[pixel];
        }

      if(DataModu[pixel]>0.2)
        {
	  sommeP_reel[pixel]=sommeP_reel[pixel]+plan_I13[pixel];
	  sommeP_imag[pixel]=sommeP_imag[pixel]+plan_I42[pixel];
	  supRedonQ[pixel]+=1;
        }
    }

  //cout<<"pixelQ_max="<<pixelQ_max<<endl;
  //cout<<"maxQholo="<<DataModu[pixelQ_max]<<endl;
  moyQholo=sommeQholo/N;
  cout<<"Moyen de qualite holo="<<moyQholo<<endl;

  SAV(DataModu, N, "/ramdisk/DataModu.bin", FLOAT,"a+b");
  SAV(supRedonQ, N, "/ramdisk/supRedonQ.bin", FLOAT,"a+b");

  memset(plan_reel, 0, N*8);
  memset(plan_imag, 0, N*8);
  //memset(supRedonQ, 0, N*8);

  //if(DataModu[pixelQ_max]>0.8)
  if(moyQholo>0.2)
    {
      for(int pixel=0;pixel<N;pixel++)
        {
	  if(supRedonQ[pixel]==0)
            {
	      supRedonQ[pixel]=1;
            }
	  plan_reel[pixel]=sommeP_reel[pixel]*(masque[pixel])/supRedonQ[pixel];
	  plan_imag[pixel]=sommeP_imag[pixel]*(masque[pixel])/supRedonQ[pixel];

	  //ampli_ref[pixel]= sqrt((double)mod_ref[pixel]);
	  //plan_reel[pixel]=sommeP_reel[pixel]*(masque[pixel])/supRedonQ[pixel]/ampli_ref[pixel];
	  //plan_imag[pixel]=sommeP_imag[pixel]*(masque[pixel])/supRedonQ[pixel]/ampli_ref[pixel];
        }
    } //fin de test de validite du nom de fichier
  else
    {
      printf("qualité d'hologramme %i insuffisante\n",cpt_angle);
    }
  //SAV(plan_reel, N, "/ramdisk/plan_reel.bin", FLOAT,"a+b");
  delete[] DataModu,sommeP_reel,sommeP_imag,plan_I42,plan_I13;
  //delete[] supRedonQ;
}

void Plan_ds_Vol(double *Vol3DRe, double *Vol3DIm, double *plan2DRe, double *plan2DIm, Var3D dimVol, Var2D dimPlan, int z3Di)
{
  int x3Di=0, y3Di=0;
  int taillePlan=dimPlan.x*dimPlan.y;
  int planVol3D=dimVol.x*dimVol.y;
  int ecartX=(dimVol.x-dimPlan.x)/2, ecartY=(dimVol.y-dimPlan.y)/2;
  int altitude_k=(z3Di+round(dimVol.z/2))*planVol3D;
  int cpt2D=0, k=0;

  for (int yi = 0; yi < dimPlan.y; yi++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
    for (int xi = 0; xi < dimPlan.x; xi++) { //on balaye l'image 2D en y, centre au milieu
      cpt2D=yi*dimPlan.x+xi;//calcul du cpt du tableau 1D de l'image 2D
      x3Di=xi+ecartX;
      y3Di=yi+ecartY;
      k=altitude_k+y3Di*dimVol.x+x3Di;
      Vol3DRe[k]=plan2DRe[cpt2D];
      Vol3DIm[k]=plan2DIm[cpt2D];
    }
  }
}
void Som_Plan_ds_Vol(double *Vol3DRe, double *Vol3DIm, double *plan2DRe, double *plan2DIm, Var3D dimVol, Var2D dimPlan, int z3Di)
{
  int x3Di=0, y3Di=0;
  int taillePlan=dimPlan.x*dimPlan.y;
  int planVol3D=dimVol.x*dimVol.y;
  int ecartX=(dimVol.x-dimPlan.x)/2, ecartY=(dimVol.y-dimPlan.y)/2;
  int altitude_k=(z3Di+round(dimVol.z/2))*planVol3D;
  int cpt2D=0, k=0;

  for (int yi = 0; yi < dimPlan.y; yi++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
    for (int xi = 0; xi < dimPlan.x; xi++) { //on balaye l'image 2D en y, centre au milieu
      cpt2D=yi*dimPlan.x+xi;//calcul du cpt du tableau 1D de l'image 2D
      x3Di=xi+ecartX;
      y3Di=yi+ecartY;
      k=altitude_k+y3Di*dimVol.x+x3Di;
      Vol3DRe[k]=plan2DRe[cpt2D]+Vol3DRe[k];
      Vol3DIm[k]=plan2DIm[cpt2D]+Vol3DIm[k];
    }
  }
}
void  retroPropagSA(float deltaZ, double *fft_shift_normRe, double *fft_shift_normIm, double * spectrePropag_normRe, double * spectrePropag_normIm,
                    Var2D dimSpctHolo, double rayon, Var2D CoordSpecI)
{
  deltaZ=deltaZ/20;
  Var2D KMAX={round(dimSpctHolo.x/2),round(dimSpctHolo.y/2)};
  const int NbPixHolo=dimSpctHolo.x*dimSpctHolo.y;
  double kz[dimSpctHolo.x*dimSpctHolo.y];
  nbCplx rephase[dimSpctHolo.x*dimSpctHolo.y];
  nbCplx spectrePropag[dimSpctHolo.x*dimSpctHolo.y];
  Var2D CoordSpecC={CoordSpecI.x-KMAX.x,CoordSpecI.y-KMAX.y};

  // nbCplx planObjet_I[dimSpctHolo.x*dimSpctHolo.y];


  ///--------création de variable pour éviter N calculs dans la boucle sur le volume 3D
  double r2=rayon*rayon;
  double arg_kz_arc=0,kz_arc=0;
  //printf("round(rayon*rayon-(xm0)^2-(ym0)^2: %i\n",round(rayon*rayon-(xm0)^2-(ym0)^2));

  double kz0 = sqrt(rayon*rayon-CoordSpecC.x*CoordSpecC.x-CoordSpecC.y*CoordSpecC.y);
  //cout<<"kz0="<<kz0<<endl;
  int NXMAX_CARRE=KMAX.x*KMAX.x;
  ///CALCULER LES KZ-----------------------------------
  for (int ky = -KMAX.y; ky < KMAX.y; ky++) { //on balaye l'image 2D en x , origine (0,0) de l'image au milieu
    int ky_carre=ky*ky;
    for (int kx = -KMAX.x; kx < KMAX.x; kx++) { //on balaye l'image 2D en y, centre au milieu
      int cpt=(ky+KMAX.y)*dimSpctHolo.x+kx+KMAX.x;//kyi*dim+kxi ...calcul du cpt du tableau 1D de l'image 2D
      if(kx*kx+ky_carre<NXMAX_CARRE)//ne pas depasser l'ouverture numérique pour 1 hologramme
	{
	  double kz_carre=r2-kx*kx-ky_carre; //altitude au carré des données
	  kz[cpt]=sqrt(kz_carre)-kz0;
	  //kz[cpt]=sqrt(kz_carre);
	}
      else  {

	kz[cpt]=0;
      }
    } //fin for y
  }
  //SAV(kz, NbPixHolo, "/home/hui/maniptomo/IDP/SA/kz.bin", FLOAT,"a+b");
  ///---------------------CALCULER LE SPECTRE REPROPAGE-----------------------------------
  double energie=0;
  for(int pix=0;pix<NbPixHolo;pix++)
    {
      rephase[pix].Re=cos(kz[pix]*deltaZ);
      rephase[pix].Im=sin(kz[pix]*deltaZ);

      spectrePropag[pix].Re=fft_shift_normRe[pix]*rephase[pix].Re-fft_shift_normIm[pix]*rephase[pix].Im;
      spectrePropag[pix].Im=fft_shift_normRe[pix]*rephase[pix].Im+fft_shift_normIm[pix]*rephase[pix].Re;
      // energie=sqrt(pow(spectrePropag[pix].Re,2)+pow(spectrePropag[pix].Im,2))+energie;
    }
  // cout<<"energie moyenne par pixel="<<energie/taillePlan<<endl;
  //SAV_Re(rephase, 4*KMAX.x*KMAX.x, "/home/mat/tomo_test/rephase.bin", FLOAT,"a+b");
  //SAV_Re(spectrePropag, NbPixHolo, "/home/mat/tomo_test/spectrePropag.bin", FLOAT,"a+b");
  //decal2DCplxGen(spectrePropag,spectrePropag_I,dimSpctHolo,KMAX);
  //SAV_Re(spectrePropag, 4*KMAX.x*KMAX.x, "/home/mat/tomo_test/spectrePropag.bin", FLOAT,"a+b");

  ///--------Renormalisation de la phase du fond après repropagation
  ///le spéculaire étant centré, on indique les coordonénes du milieu de l'image

  //cout<<"kmax="<<KMAX.x<<endl;
  int cptSpec=(dimSpctHolo.x)*CoordSpecI.y+CoordSpecI.x;
  //cout<<"cptSpec="<<cptSpec<<endl;
  //cout<<"valeur="<<spectrePropag[cptSpec].Re<<endl;

  double max_module=pow(spectrePropag[cptSpec].Re,2)+pow(spectrePropag[cptSpec].Im,2),
    max_part_reel=spectrePropag[cptSpec].Re,max_part_imag=spectrePropag[cptSpec].Im;
  for(int cpt=0; cpt<NbPixHolo; cpt++) {
    spectrePropag_normRe[cpt]=(spectrePropag[cpt].Re*max_part_reel+spectrePropag[cpt].Im*max_part_imag)/max_module;
    spectrePropag_normIm[cpt]=(spectrePropag[cpt].Im*max_part_reel-spectrePropag[cpt].Re*max_part_imag)/max_module;
    //spectrePropag_norm[cpt].Re=spectrePropag[cpt].Re;
    //spectrePropag_norm[cpt].Im=spectrePropag[cpt].Im;
  }

  // SAV_Re(spectrePropag_norm, NbPixHolo, "/home/mat/tomo_test/spectrePropag_Norm.bin", FLOAT,"a+b");
}

void Grad(double* entreeRe, double* entreeIm, double* Gradient, int taille, int NbPlan, int cptPlan)
{
  ///calcul amplitude
  double *ampli=new double[taille*taille];
  double *ampliGrad=new double[taille*taille];
  for(int pixel=0; pixel<taille*taille; pixel++)
    {
      ampli[pixel]=sqrt(pow(entreeRe[pixel],2)+pow(entreeIm[pixel],2));
    }
  //SAV(ampli, taille*taille, "/home/hui/maniptomo/IDP/SA/ampli3D.bin", FLOAT,"a+b");

  double *phase=new double[taille*taille];
  double *phaseGrad=new double[taille*taille];

  ///calcul phase de -PI à PI
  for(int pixel=0;pixel<taille*taille;pixel++)
    {
      double denominPhase=entreeRe[pixel];
      double numeraPhase=entreeIm[pixel];

      if(denominPhase>0)
        {
	  double racine=sqrt(pow(numeraPhase,2)+pow(denominPhase,2));
	  phase[pixel]=asin(numeraPhase/racine);
        }
      else
        {
	  double racine=sqrt(pow(numeraPhase,2)+pow(denominPhase,2));
	  phase[pixel]=PI-asin(numeraPhase/racine);
        }
    }
  //SAV(phase, taille*taille, "/home/hui/maniptomo/IDP/champ2d/phase3D.bin", FLOAT,"a+b");

  ///Copier double vers matrice de opencv///
  Mat img(taille,taille,CV_64F);
  for(int y=0; y<taille; y++)
    {
      for(int x=0; x<taille; x++)
	{
	  int cpt=y*taille+x;
	  img.at<double>(y,x)=ampli[cpt];
	}
    }

  /// sobel
  int scale = 1;
  int delta = 0;
  int ddepth = CV_64F;

  Mat grad_x, grad_y, imgResult;

  /// Gradient X
  //Scharr( img1[0], grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
  Sobel( img, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
  /// Gradient Y
  //Scharr( img1[0], grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
  Sobel( img, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
  /// Total Gradient (approximate)
  addWeighted( grad_x, 0.5, grad_y, 0.5, 0, imgResult );

  //    /// laplacian
  //    Mat imgResult;
  //    int kernel_size = 3;
  //    int scale = 1;
  //    int delta = 0;
  //    int ddepth = CV_64F;
  //    Laplacian(img, imgResult, ddepth, kernel_size, scale, delta, BORDER_DEFAULT );

  ///Copier matrice de opencv vers double///
  for(int y=0; y<imgResult.rows; y++)
    {
      for(int x=0; x<imgResult.cols; x++)
	{
	  int cpt=y*imgResult.cols+x;
	  ampliGrad[cpt]=imgResult.at<double>(y,x);
	}
    }
  //SAV(ampliGrad, taille*taille, "/home/hui/maniptomo/IDP/SA/ampli3D_Grad.bin", FLOAT,"a+b");

  double *ampli1= new double[taille*taille];
  double *ampli2= new double[taille*taille];
  Var2D decal1={0,-1};
  Var2D decal2={-1,0};
  Var2D dimFin={taille, taille};

  circshift2(ampliGrad,ampli1,dimFin,decal1);
  circshift2(ampliGrad,ampli2,dimFin,decal2);

  double *Grad=new double[taille*taille];
  double sommegrad=0;
  for(int cpt=0;cpt<taille*taille;cpt++)
    {
      Grad[cpt]=sqrt(pow(ampliGrad[cpt]-ampli1[cpt],2)+pow(ampliGrad[cpt]-ampli2[cpt],2));
      sommegrad=sommegrad+Grad[cpt];
    }
  //SAV(Grad, taille*taille, "/home/hui/maniptomo/IDP/champ_2d/Grad.bin", FLOAT,"a+b");

  double *sommeGrad=new double[NbPlan];
  for(int cpt=0;cpt<NbPlan;cpt++)
    {
      sommeGrad[cpt]=sommegrad;
    }
  Gradient[cptPlan]=sommeGrad[cptPlan];
  //cout<<"cptPlan= "<<cptPlan+1<<" sommegradient="<<sommegrad<<endl;
}




void AutofocusGrad(double *fft_shift_normRe, double *fft_shift_normIm, double *sortieRe, double *sortieIm, Var2D dimInit, double rayon, Var2D CoordSpecI, int cpt_angle)
{

  int zInf=-20, zSup=20, NbPlanZ=zSup-zInf;
  int xm0=(CoordSpecI.x-dimInit.x/2), ym0=(CoordSpecI.y-dimInit.x/2);
  Var3D dimVol={dimInit.x,dimInit.y,NbPlanZ};
  Var2D plandecal={dimInit.x/2, dimInit.x/2};

  double *spectrePropag_normRe=new double[dimInit.x*dimInit.x];
  double *spectrePropag_normIm=new double[dimInit.x*dimInit.x];
  double *spectrePropag_normRe_shift=new double[dimInit.x*dimInit.x];
  double *spectrePropag_normIm_shift=new double[dimInit.x*dimInit.x];

  double *planObjetRe=new double[dimInit.x*dimInit.x];
  double *planObjetIm=new double[dimInit.x*dimInit.x];
  double *planObjetRe_shift=new double[dimInit.x*dimInit.x];
  double *planObjetIm_shift=new double[dimInit.x*dimInit.x];

  double *objet3DRe=new double[dimInit.x*dimInit.x*NbPlanZ];
  double *objet3DIm=new double[dimInit.x*dimInit.x*NbPlanZ];

  double *objet3DRe_refocal=new double[dimInit.x*dimInit.x];
  double *objet3DIm_refocal=new double[dimInit.x*dimInit.x];

  double *spectreRe_refocal=new double[dimInit.x*dimInit.x];
  double *spectreIm_refocal=new double[dimInit.x*dimInit.x];

  double *Gradient= new double[NbPlanZ];

  for(int deltaZ=zInf;deltaZ<zSup;deltaZ++)
    {
      retroPropagSA(deltaZ, fft_shift_normRe, fft_shift_normIm, spectrePropag_normRe, spectrePropag_normIm, dimInit, rayon, CoordSpecI );

      Var2D decal={-xm0+dimInit.x/2,-ym0+dimInit.x/2};  ///Eliminer les franges d'inclinaision
      circshift2(spectrePropag_normRe,spectrePropag_normRe_shift,dimInit,decal);
      circshift2(spectrePropag_normIm,spectrePropag_normIm_shift,dimInit,decal);

      TF2D_INV(spectrePropag_normRe_shift,spectrePropag_normIm_shift, planObjetRe, planObjetIm, dimInit.x, dimInit.y);

      circshift2(planObjetRe,planObjetRe_shift,dimInit,plandecal);
      circshift2(planObjetIm,planObjetIm_shift,dimInit,plandecal);

      Plan_ds_Vol(objet3DRe, objet3DIm, planObjetRe_shift, planObjetIm_shift, dimVol, dimInit, deltaZ);
      Grad(planObjetRe_shift, planObjetIm_shift, Gradient, dimInit.x, NbPlanZ, deltaZ+zSup);
    }
  delete[] planObjetRe, planObjetIm, planObjetRe_shift, planObjetIm_shift;

  ///############################################
  int cpt_min=0;

  cout << endl;
  for(int cpt=0;cpt<NbPlanZ;cpt++)
    {
      cout << " | " << Gradient[ cpt ];
      if(Gradient[cpt] < Gradient[cpt_min])
	{
	  Gradient[cpt_min]=Gradient[cpt];
	  cpt_min=cpt;
	}
    }
  // triche
  cpt_min = NbPlanZ - 1;

  cout<<"cpt_min="<<cpt_min+1<<","<<"Gradient="<<Gradient[cpt_min]<<endl;
  //cout<<"cpt_max="<<cpt_max+1<<","<<"Gradient="<<Gradient[cpt_max]<<endl;
  delete[] Gradient;

  double *objet3DReTmp=new double[dimInit.x*dimInit.x*NbPlanZ];
  double *objet3DImTmp=new double[dimInit.x*dimInit.x*NbPlanZ];
  //memset(objet3DReTmp, 0, dimInit.x*dimInit.x*8);
  //memset(objet3DImTmp, 0, dimInit.x*dimInit.x*8);
  for(int cpt=0; cpt<dimInit.x*dimInit.x*NbPlanZ; cpt++)
    {
      objet3DReTmp[cpt]=objet3DRe[cpt];
      objet3DImTmp[cpt]=objet3DIm[cpt];
    }
  SAV(objet3DReTmp, dimInit.x*dimInit.x*NbPlanZ, "/ramdisk/Objet3DPlanRepro.bin", FLOAT,"a+b");

  FILE* fichier_objet3DRe = NULL;
  FILE* fichier_objet3DIm = NULL;
  fichier_objet3DRe = fopen("/ramdisk/objet3DRe_refocal.bin", "wb");
  fichier_objet3DIm = fopen("/ramdisk/objet3DIm_refocal.bin", "wb");
  MSG_ASSERT( fichier_objet3DRe , "file not found");
  MSG_ASSERT( fichier_objet3DIm , "file not found");
  int cpt=0;
  float precision;
  for(int z=cpt_min; z < cpt_min+1; z++)
    {
      for(int y=0; y < dimInit.x; y++)
	{
	  for(int x=0; x < dimInit.x; x++)
	    {
	      cpt= x + y * dimInit.x + z * dimInit.x * dimInit.x;
	      precision=objet3DReTmp[cpt];
	      fwrite(&precision,sizeof(precision),1,fichier_objet3DRe);
	      precision=objet3DImTmp[cpt];
	      fwrite(&precision,sizeof(precision),1,fichier_objet3DIm);
	    }
	}
    }
  fclose(fichier_objet3DRe);
  fclose(fichier_objet3DIm);
  //printf("ecriture de objet3DRe_refocal.bin et objet3DIm_refocal.bin : OK \n");
  delete[] objet3DReTmp, objet3DImTmp, objet3DRe, objet3DIm;
  delete[] spectrePropag_normRe, spectrePropag_normIm;

  Mat *Re= new Mat[NbPlanZ];
  Mat *Im= new Mat[NbPlanZ];
  read_dataM("/ramdisk/objet3DRe_refocal.bin", Re, dimInit.x, NbPlanZ);
  read_dataM("/ramdisk/objet3DIm_refocal.bin", Im, dimInit.x, NbPlanZ);

  ///Copier matrice de opencv vers double///
  for(int y=0; y<dimInit.x; y++)
    {
      for(int x=0; x<dimInit.x; x++)
	{
	  int cpt=y*dimInit.x+x;
	  objet3DRe_refocal[cpt]=Re[0].at<double>(y,x);
	  objet3DIm_refocal[cpt]=Im[0].at<double>(y,x);
	}
    }
  delete[] Re, Im;

  SAV(objet3DRe_refocal, dimInit.x*dimInit.y, "/ramdisk/objet3DPlanRecal.bin", FLOAT,"a+b");

  ////VECTRA: l'objet est pas circshifté en sortie
  //// conclusion: il faut le circshifter là, et donc dans l'espace image, avant le repassage freq

  double *refocus_Re = new double[dimInit.x * dimInit.x];
  double *refocus_Im = new double[dimInit.x * dimInit.x];

  ////ajout vectra
  Var2D decalC={dimInit.x/2, dimInit.x/2};
  circshift2(objet3DRe_refocal, refocus_Re, dimInit, decalC);
  circshift2(objet3DIm_refocal, refocus_Im, dimInit, decalC);
  

  TF2D(refocus_Re, refocus_Im, spectreRe_refocal, spectreIm_refocal, dimInit.x, dimInit.y);
  //   TF2D(objet3DRe_refocal,objet3DIm_refocal, spectreRe_refocal, spectreIm_refocal, dimInit.x, dimInit.y);
  delete[] objet3DRe_refocal, objet3DIm_refocal;
  
  ////ajout vectra
  delete[] refocus_Re;
  delete[] refocus_Im;


  Var2D decalI={xm0-dimInit.x/2,ym0-dimInit.x/2};
  circshift2(spectreRe_refocal, sortieRe, dimInit, decalI); ///Ajouter les franges d'inclinaision
  circshift2(spectreIm_refocal, sortieIm, dimInit, decalI);
  
    
  delete[] spectreRe_refocal, spectreIm_refocal;

}
void AutofocusMod(double *fft_shift_normRe, double *fft_shift_normIm, double *sortieRe, double *sortieIm, Var2D dimInit, double rayon, Var2D CoordSpecI, int cpt_angle)
{

  int zInf=-5, zSup=5, NbPlanZ=zSup-zInf;
  int xm0=(CoordSpecI.x-dimInit.x/2), ym0=(CoordSpecI.y-dimInit.x/2);
  Var3D dimVol={dimInit.x,dimInit.y,NbPlanZ};
  Var2D plandecal={dimInit.x/2, dimInit.x/2};

  double *spectrePropag_normRe=new double[dimInit.x*dimInit.x];
  double *spectrePropag_normIm=new double[dimInit.x*dimInit.x];
  double *spectrePropag_normRe_shift=new double[dimInit.x*dimInit.x];
  double *spectrePropag_normIm_shift=new double[dimInit.x*dimInit.x];

  double *planObjetRe=new double[dimInit.x*dimInit.x];
  double *planObjetIm=new double[dimInit.x*dimInit.x];
  double *planObjetRe_shift=new double[dimInit.x*dimInit.x];
  double *planObjetIm_shift=new double[dimInit.x*dimInit.x];

  double *objet3DRe=new double[dimInit.x*dimInit.x*NbPlanZ];
  double *objet3DIm=new double[dimInit.x*dimInit.x*NbPlanZ];

  double *objet3DRe_refocal=new double[dimInit.x*dimInit.x];
  double *objet3DIm_refocal=new double[dimInit.x*dimInit.x];

  double *spectreRe_refocal=new double[dimInit.x*dimInit.x];
  double *spectreIm_refocal=new double[dimInit.x*dimInit.x];

  double *planObjetMod_shift= new double[dimInit.x*dimInit.x];
  double *sommePlanMod= new double[NbPlanZ];

  for(int deltaZ=zInf;deltaZ<zSup;deltaZ++)
    {
      retroPropagSA(deltaZ, fft_shift_normRe, fft_shift_normIm, spectrePropag_normRe, spectrePropag_normIm, dimInit, rayon, CoordSpecI );

      Var2D decal={-xm0+dimInit.x/2,-ym0+dimInit.x/2};  ///Eliminer les franges d'inclinaision
      circshift2(spectrePropag_normRe,spectrePropag_normRe_shift,dimInit,decal);
      circshift2(spectrePropag_normIm,spectrePropag_normIm_shift,dimInit,decal);

      TF2D_INV(spectrePropag_normRe_shift,spectrePropag_normIm_shift, planObjetRe, planObjetIm, dimInit.x, dimInit.y);

      circshift2(planObjetRe,planObjetRe_shift,dimInit,plandecal);
      circshift2(planObjetIm,planObjetIm_shift,dimInit,plandecal);

      Plan_ds_Vol(objet3DRe, objet3DIm, planObjetRe_shift, planObjetIm_shift, dimVol, dimInit, deltaZ);
      for(int pixel=0; pixel<dimInit.x*dimInit.y; pixel++)
	{
	  planObjetMod_shift[pixel]=sqrt(pow(planObjetRe_shift[pixel],2)+pow(planObjetIm_shift[pixel],2));
	  sommePlanMod[deltaZ+zSup]=sommePlanMod[deltaZ+zSup]+planObjetMod_shift[pixel];
	}
    }
  delete[] planObjetRe, planObjetIm, planObjetRe_shift, planObjetIm_shift;
  delete[] planObjetMod_shift;

  ///############################################
  int cpt_min=0;
  if(cpt_angle==0)
    {
      for(int cpt=0;cpt<NbPlanZ;cpt++)
	{
	  if(sommePlanMod[cpt] < sommePlanMod[cpt_min])
	    {
	      sommePlanMod[cpt_min]=sommePlanMod[cpt];
	      cpt_min=cpt;
	    }
	}
    }
  else
    {
      for(int cpt=0;cpt<NbPlanZ;cpt++)
	{
	  if(sommePlanMod[cpt] < sommePlanMod[cpt_min])
	    {
	      sommePlanMod[cpt_min]=sommePlanMod[cpt];
	      cpt_min=cpt;
	    }
	}
    }
  cout<<"cpt_min="<<cpt_min+1<<","<<"sommePlanMod="<<sommePlanMod[cpt_min]<<endl;
  //cout<<"cpt_max="<<cpt_max+1<<","<<"sommePlanMod="<<sommePlanMod[cpt_max]<<endl;
  delete[] sommePlanMod;

  double *objet3DReTmp=new double[dimInit.x*dimInit.x*NbPlanZ];
  double *objet3DImTmp=new double[dimInit.x*dimInit.x*NbPlanZ];
  //memset(objet3DReTmp, 0, dimInit.x*dimInit.x*8);
  //memset(objet3DImTmp, 0, dimInit.x*dimInit.x*8);
  for(int cpt=0; cpt<dimInit.x*dimInit.x*NbPlanZ; cpt++)
    {
      objet3DReTmp[cpt]=objet3DRe[cpt];
      objet3DImTmp[cpt]=objet3DIm[cpt];
    }
  SAV(objet3DReTmp, dimInit.x*dimInit.x*NbPlanZ, "/ramdisk/Objet3DPlanRepro.bin", FLOAT,"a+b");

  FILE* fichier_objet3DRe = NULL;
  FILE* fichier_objet3DIm = NULL;
  fichier_objet3DRe = fopen("/ramdisk/objet3DRe_refocal.bin", "wb");
  fichier_objet3DIm = fopen("/ramdisk/objet3DIm_refocal.bin", "wb");
  MSG_ASSERT( fichier_objet3DRe , "file not found");
  MSG_ASSERT( fichier_objet3DIm , "file not found");
  int cpt=0;
  float precision;
  for(int z=cpt_min; z < cpt_min+1; z++)
    {
      for(int y=0; y < dimInit.x; y++)
	{
	  for(int x=0; x < dimInit.x; x++)
	    {
	      cpt= x + y * dimInit.x + z * dimInit.x * dimInit.x;
	      precision=objet3DReTmp[cpt];
	      fwrite(&precision,sizeof(precision),1,fichier_objet3DRe);
	      precision=objet3DImTmp[cpt];
	      fwrite(&precision,sizeof(precision),1,fichier_objet3DIm);
	    }
	}
    }
  fclose(fichier_objet3DRe);
  fclose(fichier_objet3DIm);
  //printf("ecriture de objet3DRe_refocal.bin et objet3DIm_refocal.bin : OK \n");
  delete[] objet3DReTmp, objet3DImTmp, objet3DRe, objet3DIm;
  delete[] spectrePropag_normRe, spectrePropag_normIm;

  Mat *Re= new Mat[NbPlanZ];
  Mat *Im= new Mat[NbPlanZ];
  read_dataM("/ramdisk/objet3DRe_refocal.bin", Re, dimInit.x, NbPlanZ);
  read_dataM("/ramdisk/objet3DIm_refocal.bin", Im, dimInit.x, NbPlanZ);

  ///Copier matrice de opencv vers double///
  for(int y=0; y<dimInit.x; y++)
    {
      for(int x=0; x<dimInit.x; x++)
	{
	  int cpt=y*dimInit.x+x;
	  objet3DRe_refocal[cpt]=Re[0].at<double>(y,x);
	  objet3DIm_refocal[cpt]=Im[0].at<double>(y,x);
	}
    }
  delete[] Re, Im;

  SAV(objet3DRe_refocal, dimInit.x*dimInit.y, "/ramdisk/objet3DPlanRecal.bin", FLOAT,"a+b");
  TF2D(objet3DRe_refocal,objet3DIm_refocal, spectreRe_refocal, spectreIm_refocal, dimInit.x, dimInit.y);
  delete[] objet3DRe_refocal, objet3DIm_refocal;

  Var2D decalI={xm0-dimInit.x/2,ym0-dimInit.x/2};
  circshift2(spectreRe_refocal, sortieRe, dimInit, decalI); ///Ajouter les franges d'inclinaision
  circshift2(spectreIm_refocal, sortieIm, dimInit, decalI);
  delete[] spectreRe_refocal, spectreIm_refocal;

}

double* alloPlan(int dimx, int dimy)
{
  return new double[dimx*dimy];
}

double* specSynMax(double *planSynRe, double *planSynIm, int NbPlanZ, Var2D dimFin)
{
  Mat *Re= new Mat[NbPlanZ];
  Mat *Im= new Mat[NbPlanZ];
  Mat *Amp= new Mat[NbPlanZ];
  read_dataM("/ramdisk/spectreI_reel.bin", Re, dimFin.x, NbPlanZ);
  read_dataM("/ramdisk/champ_2d/spectreI_imag.bin", Im, dimFin.x, NbPlanZ);
  read_dataM("/ramdisk/spectreI_ampli.bin", Amp, dimFin.x, NbPlanZ);

  double *spectreSynRe= new double[dimFin.x*dimFin.y];
  double *spectreSynIm= new double[dimFin.x*dimFin.y];
  double *spectreSynAmp= new double[dimFin.x*dimFin.y];

  int index_max=0;
  for(int y=0; y<dimFin.y; y++)
    {
      for(int x=0; x<dimFin.x; x++)
        {
	  for (int index = 0; index < NbPlanZ; index++)
            {
	      if(Re[index].at<double>(y,x) > Re[index_max].at<double>(y,x))
                {
		  index_max = index;
                }
	      else
                {
		  index_max = index_max;
                }
            }
	  int cpt=y*dimFin.y+x;
	  spectreSynRe[cpt]=Re[index_max].at<double>(y,x);
	  //spectreSynIm[cpt]=Im[index_max].at<double>(y,x);
	  //spectreSynAmp[cpt]=Amp[index_max].at<double>(y,x);
        }
    }

  int index_max1=0;
  for(int y=0; y<dimFin.y; y++)
    {
      for(int x=0; x<dimFin.x; x++)
        {
	  for (int index = 0; index < NbPlanZ; index++)
            {
	      if(Im[index].at<double>(y,x) > Im[index_max1].at<double>(y,x))
                {
		  index_max1 = index;
                }
	      else
                {
		  index_max1 = index_max1;
                }
            }
	  int cpt=y*dimFin.y+x;
	  //spectreSynRe[cpt]=Re[index_max].at<double>(y,x);
	  spectreSynIm[cpt]=Im[index_max1].at<double>(y,x);
	  //spectreSynAmp[cpt]=Amp[index_max].at<double>(y,x);
        }
    }

  //SAV(spectreSynRe, dimFin.x*dimFin.y, "/home/hui/maniptomo/IDP/champ_2d/spectre_tous.bin", FLOAT,"a+b");
  double *spectreSynRe_shift= new double[dimFin.x*dimFin.y];
  double *spectreSynIm_shift= new double[dimFin.x*dimFin.y];
  Var2D fftdecal={dimFin.x/2, dimFin.y/2};
  circshift2(spectreSynRe,spectreSynRe_shift,dimFin,fftdecal);
  circshift2(spectreSynIm,spectreSynIm_shift,dimFin,fftdecal);
  SAV(spectreSynRe_shift, dimFin.x*dimFin.y, "/ramdisk/spectreSynMax.bin", FLOAT,"a+b");

  TF2D_INV(spectreSynRe, spectreSynIm, planSynRe, planSynIm, dimFin.x, dimFin.y);
  SAV(planSynRe, dimFin.x*dimFin.y, "/ramdisk/PlanSynMaxRe.bin", FLOAT,"a+b");
  SAV(planSynIm, dimFin.x*dimFin.y, "/ramdisk/PlanSynMaxIm.bin", FLOAT,"a+b");

  delete[] Re, Im;
  delete[] spectreSynRe, spectreSynIm, spectreSynRe_shift, spectreSynIm_shift;
}

//void Autofocus(double *fft_shift_normRe, double *fft_shift_normIm, double *sortieRe, double *sortieIm, Var2D dimInit, double rayon, Var2D CoordSpecI)
// {
//
//     int zInf=-5, zSup=5, NbPlanZ=10;
//     int xm0=(CoordSpecI.x-dimInit.x/2), ym0=(CoordSpecI.y-dimInit.x/2);
//     Var3D dimVol={dimInit.x,dimInit.y,NbPlanZ};
//     Var2D plandecal={dimInit.x/2, dimInit.x/2};
//
//     double *spectrePropag_normRe=new double[dimInit.x*dimInit.x];
//     double *spectrePropag_normIm=new double[dimInit.x*dimInit.x];
//     double *spectrePropag_normRe_shift=new double[dimInit.x*dimInit.x];
//     double *spectrePropag_normIm_shift=new double[dimInit.x*dimInit.x];
//
//     double *planObjetRe=new double[dimInit.x*dimInit.x];
//     double *planObjetIm=new double[dimInit.x*dimInit.x];
//     double *planObjetRe_shift=new double[dimInit.x*dimInit.x];
//     double *planObjetIm_shift=new double[dimInit.x*dimInit.x];
//
//     double** objet3DRe;
//     objet3DRe=new double*[NbPlanZ];
//     double** objet3DIm;
//     objet3DIm=new double*[NbPlanZ];
//     for(int cpt=0; cpt< NbPlanZ; cpt++)
//     {
//         objet3DRe[cpt] = alloPlan(dimInit.x, dimInit.y);
//         objet3DIm[cpt] = alloPlan(dimInit.x, dimInit.y);
//     }
//
//
//     double *Gradient= new double[NbPlanZ];
//
//     for(int deltaZ=zInf;deltaZ<zSup;deltaZ++)
//     {
//         retroPropagSA(deltaZ, fft_shift_normRe, fft_shift_normIm, spectrePropag_normRe, spectrePropag_normIm, dimInit, rayon, CoordSpecI );
//
//         Var2D decal={-xm0+dimInit.x/2,-ym0+dimInit.x/2};
//         circshift2(spectrePropag_normRe,spectrePropag_normRe_shift,dimInit,decal);
//         circshift2(spectrePropag_normIm,spectrePropag_normIm_shift,dimInit,decal);
//
//         TF2D_INV(spectrePropag_normRe_shift,spectrePropag_normIm_shift, planObjetRe, planObjetIm, dimInit.x, dimInit.y);
//
//         circshift2(planObjetRe,planObjetRe_shift,dimInit,plandecal);
//         circshift2(planObjetIm,planObjetIm_shift,dimInit,plandecal);
//
//         Plan_ds_Vol(objet3DRe[deltaZ], objet3DIm[deltaZ], planObjetRe_shift, planObjetIm_shift, dimVol, dimInit, deltaZ);
//         Grad(planObjetRe_shift, planObjetIm_shift, Gradient, dimInit.x, NbPlanZ, deltaZ+5);
//     }
//
/////############################################
////        int cpt_max=0;
////        for(int cpt=0;cpt<NbPlanZ;cpt++)
////        {
////            if(Gradient[cpt] > Gradient[cpt_max])
////            {
////                Gradient[cpt_max]=Gradient[cpt];
////                cpt_max=cpt;
////            }
////        }
//        int cpt_min=0;
//        for(int cpt=0;cpt<NbPlanZ;cpt++)
//        {
//            if(Gradient[cpt] < Gradient[cpt_min])
//            {
//                Gradient[cpt_min]=Gradient[cpt];
//                cpt_min=cpt;
//            }
//        }
//        cout<<"cpt_min="<<cpt_min+1<<","<<"Gradient="<<Gradient[cpt_min]<<endl;
//        //cout<<"cpt_max="<<cpt_max+1<<","<<"Gradient="<<Gradient[cpt_max]<<endl;
//        delete[] Gradient;
//
//        for(int cpt=0; cpt<dimInit.x*dimInit.y; cpt++)
//        {
//            objet3DReMin[cpt]=objet3DRe[cpt_min][cpt];
//            objet3DImMin[cpt]=objet3DIm[cpt_min][cpt];
//        }
//
//        delete[] planObjetRe, planObjetIm, planObjetRe_shift, planObjetIm_shift;
//
// }

void ContrSpeckle(double* holo1, int taille, double contrSpeckle)
{
  double *Inten= new double[taille*taille];
  double sommeInten=0;
  double meanInten=0;

  for(int pixel=0; pixel< taille*taille; pixel++)
    {
      Inten[pixel]=holo1[pixel];
      sommeInten=sommeInten+Inten[pixel];
    }
  meanInten=sommeInten/(taille*taille);

  double sigma=0;
  double sommeSigma=0;
  for(int pixel=0; pixel< taille*taille; pixel++)
    {
      sommeSigma=sommeSigma+pow(Inten[pixel]-meanInten,2);
    }
  sigma=sqrt(sommeSigma/(taille*taille-1));
  //double contrSpeckle;
  contrSpeckle=sigma/meanInten;
  cout<<"Contraste de Speckle="<<contrSpeckle<<endl;
}

