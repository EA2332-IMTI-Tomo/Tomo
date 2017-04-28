/***************************************************************
 * Name:      Gui_TomoMain.h
 * Purpose:   Defines Application Frame
 * Author:    Matthieu Debailleul ()
 * Created:   2017-03-28
 * Copyright: Matthieu Debailleul ()
 * License:
 **************************************************************/
#include <iostream>

#ifndef GUI_TOMOMAIN_H
#define GUI_TOMOMAIN_H

#ifndef WX_PRECOMP
    #include <wx/wx.h>
#endif

#include <fstream>
#include <sstream>
#include <math.h>
#include<vector>

class Gui_TomoFrame: public wxFrame
{   ///La classe de la fenêtre principale est publique
    public:

        Gui_TomoFrame(wxFrame *frame, const wxString& title);
        ~Gui_TomoFrame();
         //float extract_val(std::string token,  std::string chemin_fic);
         std::string extract_string(std::string token,  std::string chemin_fic);
         void sav_val(std::string chemin_fic,std::vector<std::string> &tab_val);
         void init_tab_val(std::string chemin_fic, std::vector<std::string> &tab_val);
         void modif_tab_val(std::string token,std::string string_valeur_token,std::vector<std::string> &tab_val);

  protected:
      ///ENUM servant à identifier les boutons et menus
        enum
        {
            idMenuQuit = 1000,
            wxID_EXIT,
            idMenuAbout,
            idBoutonManip,
            idBoutonTraitement,
            idBoutonReconstruction,
            id_boolBorn,
            idBoolDeroul,
            idBoolAber,
            idBoutonOpenDir,
            idBoutonOpenDirResult,
            idBoutonRecons,
            idBoutonOpenMask,
            idBoutonRecalc,
            idBoutonSAV,
            idMenuSAV,
        };
        double Vxmax,Vymax,Vxmin,Vymin;
        int fx0,fy0,nbHolo;
        unsigned int  DIM_FINAL, NXMAX;
        bool b_BORN, b_DEROUL, b_ABER;
       // std::vector<std::string> fichier_tmp;
        std::vector<std::string> tab_val_recon;
        std::vector<std::string> tab_val_manip;
        std::vector<std::string> tab_val_gui_conf;


    // textFicManip=new wxTextCtrl(zone1,-1,"",wxPoint(160,195),wxSize(100,20));


        ///Déclaration des attributs de type Boutons et menus en protected
        std::string chemin_config_GUI,chemin_recon, chemin_config_manip,chemin_acquis,chemin_result,repertoire_config;
        wxButton* BoutonManip;
        wxButton* BoutonTraitement;
        wxButton* BoutonRecons;
        wxButton* BoutonOpenDir,*BoutonOpenDirResult,*BoutonOpenMask;
        wxBitmapButton *BoutonRecalc, *BoutonSAV;
        wxCheckBox *m_born,*m_Aber,*m_Deroul;

        wxStaticText* t;
        wxTextCtrl *editX,*editY,*editResult,*editDirAcquis,*editFicMask,*editFicManip,*editDirResultAcquis;//champ de texte editable
        wxTextCtrl *editNbHolo,*editDimFinal;
        wxTextCtrl *editCX,*editCY,*editNXMAX,*editVxmin,*editVymin,*editVxmax,*editVymax,*editDeltaVx,*editDeltaVy;
        wxStaticText *textX,*textY,*textResult,*textDirAcquis,*textFicManip;
        wxStaticText *titre_Acquis,*titre_Pretraitement, *titre_Recons,*titre_HorsAxe,*titre_Balayage;//texte des champs
        wxStaticText *textNbHolo,*textDimFinal,*textCX,*textCY,*textNXMAX,*textVxmin,*textVymin,*textVxmax,*textVymax,*textDeltaVy,*textDeltaVx,*textOffx,*textOffy;//texte des champs


  private:
        void OnClose(wxCloseEvent& event);
        void OnQuit(wxCommandEvent& event);
        void OnAbout(wxCommandEvent& event);
        void OnBoutonManip(wxCommandEvent& event);
        void OnBoutonTraitement(wxCommandEvent& event);
        void OnBoutonRecons(wxCommandEvent& event);
        void OnBoutonBorn(wxCommandEvent& event);
        void OnBoutonAber(wxCommandEvent& event);
        void OnBoutonDeroul(wxCommandEvent& event);
        void OnBoutonOpenDir(wxCommandEvent& event);
        void OnBoutonOpenDirResult(wxCommandEvent& event);
        void OnBoutonOpenMask(wxCommandEvent& event);
        //void OnBoutonRecalc(wxCommandEvent& WXUNUSED(event));
        void OnBoutonSAV(wxCommandEvent& WXUNUSED(event));

        DECLARE_EVENT_TABLE()
};


#endif // GUI_TOMOMAIN_H
