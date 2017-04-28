/***************************************************************
 * Name:      Gui_TomoMain.cpp
 * Purpose:   Code for Application Frame
 * Author:    Matthieu Debailleul ()
 * Created:   2017-03-28
 * Copyright: Matthieu Debailleul ()
 * License:
 **************************************************************/
#include "Gui_TomoApp.h"

#ifdef WX_PRECOMP
#include "wx_pch.h"
#endif

#ifdef __BORLANDC__
#pragma hdrstop
#endif //__BORLANDC__

#include "Gui_TomoMain.h"

using namespace std;

//helper functions
enum wxbuildinfoformat {
    short_f, long_f };

wxString wxbuildinfo(wxbuildinfoformat format)
{
    wxString wxbuild(wxVERSION_STRING);

    if (format == long_f )
    {
#if defined(__WXMSW__)
        wxbuild << _T("-Windows");
#elif defined(__WXMAC__)
        wxbuild << _T("-Mac");
#elif defined(__UNIX__)
        wxbuild << _T("-Linux");
#endif

#if wxUSE_UNICODE
        wxbuild << _T("-Unicode build");
#else
        wxbuild << _T("-ANSI build");
#endif // wxUSE_UNICODE
    }

    return wxbuild;
}
using namespace std;
BEGIN_EVENT_TABLE(Gui_TomoFrame, wxFrame)
    EVT_CLOSE(Gui_TomoFrame::OnClose)
    EVT_MENU(wxID_EXIT, Gui_TomoFrame::OnQuit)
    EVT_MENU(idMenuAbout, Gui_TomoFrame::OnAbout)
    EVT_BUTTON(idBoutonManip, Gui_TomoFrame::OnBoutonManip)
    EVT_BUTTON(idBoutonTraitement, Gui_TomoFrame::OnBoutonTraitement)
    EVT_BUTTON(idBoutonRecons, Gui_TomoFrame::OnBoutonRecons)
   // EVT_BUTTON(idBoutonCalcul, Gui_TomoFrame::OnBoutonCalcul)
    EVT_CHECKBOX(id_boolBorn,Gui_TomoFrame::OnBoutonBorn)
    EVT_CHECKBOX(idBoolDeroul,Gui_TomoFrame::OnBoutonDeroul)
    EVT_CHECKBOX(idBoolAber,Gui_TomoFrame::OnBoutonAber)
    EVT_BUTTON(idBoutonOpenDir, Gui_TomoFrame::OnBoutonOpenDir)
    EVT_BUTTON(idBoutonOpenDirResult, Gui_TomoFrame::OnBoutonOpenDirResult)
    EVT_BUTTON(idBoutonOpenMask, Gui_TomoFrame::OnBoutonOpenMask)
    //EVT_BUTTON(idBoutonRecalc, Gui_TomoFrame::OnBoutonRecalc)
    EVT_BUTTON(idBoutonSAV, Gui_TomoFrame::OnBoutonSAV)
    EVT_MENU(idMenuSAV, Gui_TomoFrame::OnBoutonSAV)
END_EVENT_TABLE()
///constructeur de la fenêtre, on y met toutes les init de boutons/menus (faisable aussi dans le App.main).
Gui_TomoFrame::Gui_TomoFrame(wxFrame *frame, const wxString& title)
    : wxFrame(frame, -1, title)
{
#if wxUSE_MENUS
    // create a menu bar
    wxMenuBar* mbar = new wxMenuBar();
    wxMenu* fileMenu = new wxMenu(_T(""));

    fileMenu->Append(idMenuSAV, _("&Sauver\tAlt-S"), wxT("Sauver les paramètres"));
    fileMenu->Append(wxID_EXIT, _("&Quitter\tAlt-F4"), _("Quitter l'application"));
    mbar->Append(fileMenu, _("&File"));

    wxMenu* helpMenu = new wxMenu(_T(""));
    helpMenu->Append(idMenuAbout, _("&A propos\tF1"), _("Information sur l'application"));
    mbar->Append(helpMenu, _("&Aide"));

    SetMenuBar(mbar);

#endif // wxUSE_MENUS

#if wxUSE_STATUSBAR
    // create a status bar with some information about the used wxWidgets version
    CreateStatusBar(2);
    SetStatusText(_("Tomo Mithra 0.1!"),0);
    SetStatusText(wxbuildinfo(short_f), 1);
#endif // wxUSE_STATUSBAR



///Initialiser les deux tableaux stockant les valeurs de recon.txt et config_manip.txt en mémoire

    string home=getenv("HOME");
    string fin_chemin_gui_tomo="/.config/gui_tomo.conf";
    chemin_config_GUI=getenv("HOME")+fin_chemin_gui_tomo;
    cout<<"lecture des chemins dans le fichier "<<chemin_config_GUI<<endl;

    repertoire_config=extract_string("CHEMIN_CONFIG",home+fin_chemin_gui_tomo);
    chemin_acquis=extract_string("CHEMIN_ACQUIS",home+fin_chemin_gui_tomo);
    chemin_result=extract_string("CHEMIN_RESULT",home+fin_chemin_gui_tomo);
    chemin_recon=repertoire_config+"recon.txt";
    chemin_config_manip=repertoire_config+"config_manip.txt";

    vector<string> &r_tab_val_gui_conf=tab_val_gui_conf;
    vector<string> &r_tab_val_recon=tab_val_recon;
    vector<string> &r_tab_val_manip=tab_val_manip;
    this->init_tab_val(chemin_recon,r_tab_val_recon);
    this->init_tab_val(chemin_config_manip,r_tab_val_manip);
    this->init_tab_val(home+fin_chemin_gui_tomo,r_tab_val_gui_conf);
   //BoutonTraitement= new wxButton(this, idBoutonTraitement, wxT("&Traitement"), wxPoint(50,100), wxDefaultSize, 0);

    ///text statique affiché à un certaine position
    //t = new wxStaticText(this, idBoutonTraitement, wxT("Whizzo"), wxPoint(100,150));
    //t->SetToolTip(wxT("Some incredibly compelling text"));

    //Création du wxBoxSizer horizontal (sizer principal englobant)
    wxBoxSizer *sizer_horizontal = new wxBoxSizer(wxHORIZONTAL);
    // Affectation du wxBoxSizer horizontal à la fenêtre
    SetSizer(sizer_horizontal);

            ///#################### Partie Gauche #############################
            // Création du wxPanel pour la zone de gauche
        wxBoxSizer *sizer_vertical0 = new wxBoxSizer(wxVERTICAL);
        sizer_horizontal->Add(sizer_vertical0, 1, wxALL | wxEXPAND, 1); // On l'ajoute au premier wxBoxSizer
            wxPanel *zone10 = new wxPanel(this);
            sizer_vertical0->Add(zone10, 0, wxALL | wxEXPAND, 1);
            zone10->SetBackgroundColour(*wxLIGHT_GREY);
                ///Titre
                titre_Acquis=new wxStaticText(zone10,-1,"Acquisition",wxPoint(0,0),wxSize(100,40));
                titre_Acquis->SetFont( wxFont(12, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL) );


              ///--------------------NB HOLO
            wxPanel *zone11 = new wxPanel(this);
            sizer_vertical0->Add(zone11, 0, wxALL | wxEXPAND | wxALIGN_CENTER, 1);// Ajout de cette zone au sizer vertical    //
            zone11->SetBackgroundColour(*wxLIGHT_GREY);
                textNbHolo=new wxStaticText(zone11,-1,"Nombre d'hologrammes  : ",wxPoint(10,10),wxSize(100,50));

                //édition du champ
                editNbHolo=new wxTextCtrl(zone11,-1,"",wxPoint(120,18),wxSize(40,25));
                editNbHolo->SetToolTip(wxT("Nombre d'hologrammes à acquérir"));
              //int Nb_Holo=extract_val("NB_HOLO",chemin_config_manip);
                int Nb_Holo=stof(extract_string("NB_HOLO",chemin_config_manip));//stof=string to float

                wxString string_NbHolo = wxString::Format(wxT("%i"),Nb_Holo);
                editNbHolo->SetValue(string_NbHolo);



            ///--------------------BALAYAGE------------------
            wxPanel *zone12 = new wxPanel(this);
            sizer_vertical0->Add(zone12, 0, wxALL | wxEXPAND | wxALIGN_CENTER, 1);// Ajout de cette zone au sizer vertical    //
            zone12->SetBackgroundColour(*wxLIGHT_GREY);
                    titre_Balayage=new wxStaticText(zone12,-1," Balayage",wxPoint(0,0),wxSize(100,40));
                    titre_Balayage->SetFont(wxFont(10, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD) );

                    textVxmin=new wxStaticText(zone12,-1,"Vxm ",wxPoint(10,30),wxSize(30,40));
                                editVxmin=new wxTextCtrl(zone12,-1,"",wxPoint(40,28),wxSize(50,25));
                                Vxmin=stof(extract_string("VXMIN",chemin_config_manip));
                                wxString string_Vxmin = wxString::Format(wxT("%.2f"),Vxmin);
                                editVxmin->SetValue(string_Vxmin);
                                editVxmin->SetToolTip(wxT("Tension min x"));

                    textVymin=new wxStaticText(zone12,-1,"Vym ",wxPoint(105,30),wxSize(30,40));
                                editVymin=new wxTextCtrl(zone12,-1,"",wxPoint(135,28),wxSize(50,25));
                                Vymin=stof(extract_string("VYMIN",chemin_config_manip));
                                wxString string_Vymin = wxString::Format(wxT("%.2f"),Vymin);
                                editVymin->SetValue(string_Vymin);
                                editVymin->SetToolTip(wxT("Tension min y"));

                    textVxmax=new wxStaticText(zone12,-1,"VxM ",wxPoint(10,70),wxSize(30,40));
                                editVxmax=new wxTextCtrl(zone12,-1,"",wxPoint(40,68),wxSize(50,25));
                                Vxmax=stof(extract_string("VXMAX",chemin_config_manip));
                                wxString string_Vxmax = wxString::Format(wxT("%.2f"),Vxmax);
                                editVxmax->SetValue(string_Vxmax);
                                editVxmax->SetToolTip(wxT("Tension max axe x"));

                    textVymax=new wxStaticText(zone12,-1,"VyM ",wxPoint(105,70),wxSize(30,40));
                                editVymax=new wxTextCtrl(zone12,-1,"",wxPoint(135,68),wxSize(50,25));
                                Vymax=stof(extract_string("VYMAX",chemin_config_manip));
                                wxString string_Vymax = wxString::Format(wxT("%.2f"),Vymax);
                                editVymax->SetValue(string_Vymax);
                                editVymax->SetToolTip(wxT("Tension max axe y"));



                    /* textDeltaVx=new wxStaticText(zone12,-1,wxT("vX0"),wxPoint(10,130),wxSize(25,40));
                                editDeltaVx=new wxTextCtrl(zone12,-1,"",wxPoint(40,126),wxSize(56,25),wxTE_READONLY);
                                float DeltaVx=stof(extract_string("DELTA_VX",chemin_config_manip));
                                float Rx=abs(Vxmax-Vxmin);
                                float offX=Vxmax-Rx;
                                wxString string_DeltaVx = wxString::Format(wxT("%.2f"),DeltaVx);
                                editDeltaVx->SetValue(string_DeltaVx);
                                editDeltaVx->SetToolTip(wxT("Offset tension x"));

                     textDeltaVy=new wxStaticText(zone12,-1,wxT("vY0"),wxPoint(105,130),wxSize(27,40));
                                editDeltaVy=new wxTextCtrl(zone12,-1,"",wxPoint(135,126),wxSize(54,25),wxTE_READONLY);
                                //float DeltaVy=stof(extract_string("DELTA_VY",chemin_config_manip));
                                float Ry=abs(Vymax-Vymin);
                                float offY=Vymax-Ry;
                                wxString string_DeltaVy = wxString::Format(wxT("%.2f"),DeltaVy);
                                editDeltaVy->SetValue(string_DeltaVy);
                                editDeltaVy->SetToolTip(wxT("Offset tension y"));*/


                  //  BoutonRecalc=new wxBitmapButton(zone12, idBoutonRecalc, wxBitmap(wxT("recalc.png"),
                   // wxBITMAP_TYPE_PNG), wxPoint(200, 120));


            ///--------------------ACQUISITION------------------
            wxPanel *zone13 = new wxPanel(this);
            sizer_vertical0->Add(zone13, 1, wxALL | wxEXPAND | wxALIGN_CENTER, 1);// Ajout de cette zone au sizer vertical    //
            zone13->SetBackgroundColour(*wxLIGHT_GREY);
                    /// Dialog fichier;
                    BoutonOpenDir=new wxButton(zone13, idBoutonOpenDir, wxT("&Données"), wxPoint(5,10), wxDefaultSize, 0);
                    editDirAcquis=new wxTextCtrl(zone13,-1,chemin_acquis,wxPoint(105,16), wxSize(175,26));
                    editDirAcquis->SetToolTip(wxT("Répertoire des acquisitions"));
                    // textFicManip=new wxTextCtrl(zone10,-1,"",wxPoint(160,195),wxSize(100,20));
                    BoutonManip= new wxButton(zone13, idBoutonManip, wxT("&Acquisition"), wxPoint(195,75), wxDefaultSize, 0);



    ///#################### Partie centrale : PRETRAITEMENT#############################
            //*************Création du wxBoxSizer vertical pour la partie centrale**********
                wxBoxSizer *sizer_vertical = new wxBoxSizer(wxVERTICAL);
                sizer_horizontal->Add(sizer_vertical, 1, wxALL | wxEXPAND, 1); // On l'ajoute au premier wxBoxSizer

                        //-------------zone du haut -------------------------------------------------------------------------------------
                        wxPanel *zone20 = new wxPanel(this);
                        zone20->SetBackgroundColour(*wxLIGHT_GREY);                                  //
                            sizer_vertical->Add(zone20, 0, wxALL | wxEXPAND | wxALIGN_CENTER, 1);// Ajout de cette zone au sizer vertical    //
                            titre_Pretraitement=new wxStaticText(zone20,-1,wxT(" Prétraitement"),wxPoint(0,0),wxSize(100,40));               //
                            titre_Pretraitement->SetFont(wxFont(12, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL) );      //
                        //-------------zone du milieu ------------------------------------------------------------------------------------
                        wxPanel *zone21 = new wxPanel(this, wxID_ANY, wxDefaultPosition, wxSize(-1, 100));
                        zone21->SetBackgroundColour(*wxLIGHT_GREY);

                            // Ajout de cette zone au sizer vertical
                            sizer_vertical->Add(zone21, 0, wxALL | wxEXPAND, 1);

                            m_born = new wxCheckBox(zone21, id_boolBorn, wxT("Utiliser Born"), wxPoint(20, 11));
                            m_born->SetToolTip(wxT("0=Rytov"));
                            m_born->SetValue(stof(extract_string("BORN",chemin_recon)));
                            //m_born->SetValue(true);
                            m_Deroul = new wxCheckBox(zone21, idBoolDeroul, wxT("Déroulement"), wxPoint(20, 41));
                            m_Deroul->SetValue(stof(extract_string("DEROUL",chemin_recon)));
                            m_Aber = new wxCheckBox(zone21, idBoolAber, wxT("Correction aberration"), wxPoint(20, 71));
                            m_Aber->SetValue(stof(extract_string("C_ABER",chemin_recon)));


                        //-------------zone du bas---------------------------------------------------------------------------------------//
                        wxPanel *zone22 = new wxPanel(this);
                        zone22->SetBackgroundColour(*wxLIGHT_GREY);
                            // Ajout de cette zone au sizer vertical
                            sizer_vertical->Add(zone22, 2, wxALL | wxEXPAND, 1);
                                 /// Dialog fichier;
                                BoutonOpenMask=new wxButton(zone22, idBoutonOpenMask, wxT("&Masque"), wxPoint(5,10), wxSize(70,30), 0);
                                editFicMask=new wxTextCtrl(zone22,-1,chemin_acquis+"Mask.png",wxPoint(85,6), wxSize(170,40));
                                editFicMask->SetToolTip(wxT("Masque pour correction d'aberration"));

    ///#############################Partie droite : Reconstruction#################################

                wxBoxSizer *sizer_vertical2 = new wxBoxSizer(wxVERTICAL);
                sizer_horizontal->Add(sizer_vertical2, 1, wxALL | wxEXPAND, 1);// On l'ajoute au premier wxBoxSizer

                        wxPanel *zone30= new wxPanel(this);
                        zone30->SetBackgroundColour(*wxLIGHT_GREY);
                                // Ajout de cette zone au sizer vertical
                                sizer_vertical2->Add(zone30, 0, wxALL | wxEXPAND | wxALIGN_CENTER, 1);
                                titre_Recons=new wxStaticText(zone30,-1,wxT(" Reconstruction"),wxPoint(0,0),wxSize(100,40));
                                titre_Recons->SetFont(wxFont(12, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_NORMAL) );
                        //-------------zone Hors axe -------------------------------------------------------------------------------------
                        wxPanel *zone31 = new wxPanel(this);
                        zone31->SetBackgroundColour(*wxLIGHT_GREY);
                                sizer_vertical2->Add(zone31, 0, wxALL | wxEXPAND | wxALIGN_CENTER, 1);
                                titre_HorsAxe=new wxStaticText(zone31,-1,wxT(" Hors-axe"),wxPoint(0,0),wxSize(100,40));
                                titre_HorsAxe->SetFont(wxFont(10, wxFONTFAMILY_TELETYPE, wxFONTSTYLE_NORMAL, wxFONTWEIGHT_BOLD) );

                                //CX-----
                                textCX=new wxStaticText(zone31,-1,"fx0 ",wxPoint(10,30),wxSize(26,40));
                                editCX=new wxTextCtrl(zone31,-1,"",wxPoint(36,25),wxSize(50,25));
                                editCX->SetToolTip(wxT("Centre X du spectre"));
                                int CX=stof(extract_string("CIRCLE_CX",chemin_config_manip));
                                wxString string_CX = wxString::Format(wxT("%i"),CX);
                                editCX->SetValue(string_CX);
                                //CY------
                                textCY=new wxStaticText(zone31,-1,"fy0 ",wxPoint(91,30),wxSize(28,40));
                                editCY=new wxTextCtrl(zone31,-1,"",wxPoint(115,25),wxSize(50,25));
                                int CY=stof(extract_string("CIRCLE_CY",chemin_config_manip));
                                wxString string_CY = wxString::Format(wxT("%i"),CY);
                                editCY->SetValue(string_CY);
                                editCY->SetToolTip(wxT("Centre Y du spectre"));
                                //NXMAX------
                                textNXMAX=new wxStaticText(zone31,-1,"NXMAX ",wxPoint(168,30),wxSize(60,40));
                                editNXMAX=new wxTextCtrl(zone31,-1,"",wxPoint(225,25),wxSize(50,25));
                                int NXMAX=stof(extract_string("NXMAX",chemin_config_manip));
                                wxString string_NXMAX = wxString::Format(wxT("%i"),NXMAX);
                                editNXMAX->SetValue(string_NXMAX);
                                editNXMAX->SetToolTip(wxT("Demi largeur du spectre"));

                                //DIMFINALE
                                textDimFinal=new wxStaticText(zone31,-1,"DimFinal ",wxPoint(166,65),wxSize(60,40));
                                editDimFinal=new wxTextCtrl(zone31,-1,"",wxPoint(225,60),wxSize(50,25));
                                int DimFinal=stof(extract_string("DIM_FINAL",chemin_recon));
                                wxString string_DimFinal = wxString::Format(wxT("%i"),DimFinal);
                                editDimFinal->SetValue(string_DimFinal);
                                editDimFinal->SetToolTip(wxT("Dimension du volume final"));


                        //-------------zone du bas -------------------------------------------------------------------------------------
                        wxPanel *zone32 = new wxPanel(this);
                        zone32->SetBackgroundColour(*wxLIGHT_GREY);
                                // Ajout de cette zone au sizer vertical
                                sizer_vertical2->Add(zone32, 1, wxALL | wxEXPAND | wxALIGN_CENTER, 1);


                                BoutonOpenDirResult=new wxButton(zone32, idBoutonOpenDirResult, wxT("& Résultats "), wxPoint(5,10), wxDefaultSize, 0);
                                editDirResultAcquis=new wxTextCtrl(zone32,-1,chemin_result,wxPoint(105,16),wxSize(170,20));
                                editDirResultAcquis->SetToolTip(wxT("Répertoire résultats"));
                                // textFicManip=new wxTextCtrl(zone10,-1,"",wxPoint(160,195),wxSize(100,20));
                                BoutonRecons= new wxButton(zone32, idBoutonRecons, wxT("&Reconstruction"), wxPoint(175,144), wxDefaultSize, 0);

                              //  BoutonSAV=new wxBitmapButton(zone32, idBoutonSAV, wxBitmap(wxT("save-icon.png"),
                                //wxBITMAP_TYPE_PNG), wxPoint(235, 200));
                               // BoutonSAV->SetToolTip(wxT("Sauvegarder l'ensemble des paramètres"));

}

Gui_TomoFrame::~Gui_TomoFrame()
{

}

void Gui_TomoFrame::OnClose(wxCloseEvent &event)
{
    Destroy();
}

void Gui_TomoFrame::OnQuit(wxCommandEvent &event)
{
    Destroy();
}
void Gui_TomoFrame::OnAbout(wxCommandEvent &event)
{
    wxString msg = ("Interface graphique \npour microscope tomographique \nV.0.1");
    wxMessageBox(msg, _("Bienvenue..."));

}
void Gui_TomoFrame::OnBoutonManip(wxCommandEvent &event)
{
    float Rx=abs(Vxmax-Vxmin);
    float offX=Vxmax-Rx;
    float Ry=abs(Vymax-Vymin);
    float offY=Vymax-Ry;

    wxString param1 =editDirAcquis->GetValue();
    wxString string_NbHolo=editNbHolo->GetValue();
    wxString cmd_final=wxT("maniptomo5 -ni ")+editNbHolo->GetValue()+wxT(" -voffset ")+wxString::Format(wxT("%.2f"), offX)+" "+
                        wxString::Format(wxT("%.2f"),offY)+wxT(" -vfleur ")+wxString::Format(wxT("%.2f"),Rx)+" "+wxString::Format(wxT("%.2f"),Ry);
    cout<<cmd_final<<endl;
    wxExecute(cmd_final);
}

void Gui_TomoFrame::OnBoutonTraitement(wxCommandEvent &event)
{
    wxExecute(wxT("ls -lah"));
}

void Gui_TomoFrame::OnBoutonRecons(wxCommandEvent &event)
{
   wxExecute(wxT("Tomo_pretraitement -i ")+editDirResultAcquis->GetValue());
}

void Gui_TomoFrame::OnBoutonBorn(wxCommandEvent& WXUNUSED(event))
{
  if (m_born->GetValue()==true) {
        modif_tab_val("BORN","1",tab_val_recon);
  } else {
      modif_tab_val("BORN","0",tab_val_recon);
      modif_tab_val("DEROUL", "1",tab_val_recon);
      m_Deroul->SetValue("1");
      }
}

void Gui_TomoFrame::OnBoutonDeroul(wxCommandEvent& WXUNUSED(event))
{
  if(m_Deroul->GetValue()==true) {
      modif_tab_val("DEROUL", "1",tab_val_recon);
  } else {
      modif_tab_val("DEROUL", "0",tab_val_recon);
  }
}

void Gui_TomoFrame::OnBoutonAber(wxCommandEvent& WXUNUSED(event))
{
  if (m_Aber->GetValue()==true) {
      modif_tab_val("C_ABER", "1",tab_val_recon);
  } else {
      modif_tab_val("C_ABER", "0",tab_val_recon);
  }
}

void Gui_TomoFrame::OnBoutonOpenDir(wxCommandEvent& WXUNUSED(event))
{
    wxString defaultPath = chemin_acquis;//wxT("/ramdisk/Acquis/");
    wxDirDialog* OpenDirDial = new wxDirDialog(this, _("Choisir un répertoire"),defaultPath);

    if ( OpenDirDial->ShowModal() == wxID_OK ){
            wxString path = OpenDirDial->GetPath();
            editDirAcquis->SetValue(path);
        }
}

void Gui_TomoFrame::OnBoutonOpenDirResult(wxCommandEvent& WXUNUSED(event))
{
    wxString defaultPath = chemin_acquis;
    wxDirDialog* OpenDirDial = new wxDirDialog(this, _("Choisissez un répertoire"),defaultPath);

    if ( OpenDirDial->ShowModal() == wxID_OK )
        {
            wxString path = OpenDirDial->GetPath();
            editDirResultAcquis->SetValue(path);
        }
}

void Gui_TomoFrame::OnBoutonOpenMask(wxCommandEvent& WXUNUSED(event))
{
    wxString defaultPath = chemin_acquis;
    wxFileDialog* OpenFileDial = new wxFileDialog(this, _("Choisissez un répertoire"),defaultPath);
    if ( OpenFileDial->ShowModal() == wxID_OK ){
            wxString path = OpenFileDial->GetPath();
            editFicMask->SetValue(path);
        }
}
/*
void Gui_TomoFrame::OnBoutonRecalc(wxCommandEvent& WXUNUSED(event))
{
  wxString strVymax,strVxmax,strVymin,   strVxmin;

    strVymax=editVymax->GetValue();
    strVymax.ToDouble(&Vymax);
    strVymin=editVymin->GetValue();
    strVymin.ToDouble(&Vymin);

    strVxmax=editVxmax->GetValue();
    strVxmax.ToDouble(&Vxmax);
    strVxmin=editVxmin->GetValue();
    strVxmin.ToDouble(&Vxmin);

    float VY0=(Vymax+Vymin)/2;
    float VX0=(Vxmax+Vxmin)/2;

    wxString string_VY0 = wxString::Format(wxT("%.2f"),VY0);
    editDeltaVy->SetValue(string_VY0);
    wxString string_VX0 = wxString::Format(wxT("%.2f"),VX0);
    editDeltaVx->SetValue(string_VX0);
}*/
void Gui_TomoFrame::OnBoutonSAV(wxCommandEvent& WXUNUSED(event))
{
  //tableau de config de l'interface graphique (chemins)
  modif_tab_val("CHEMIN_RESULT",editDirResultAcquis->GetValue().ToStdString(),tab_val_gui_conf);
  modif_tab_val("CHEMIN_ACQUIS",editDirAcquis->GetValue().ToStdString(),tab_val_gui_conf);
  modif_tab_val("CHEMIN_MASK",editFicMask->GetValue().ToStdString(),tab_val_gui_conf);

  //tableau de config  reconstruction
  modif_tab_val("FINAL_ANGLE",editNbHolo->GetValue().ToStdString(),tab_val_recon);
  modif_tab_val("DIM_FINAL",editDimFinal->GetValue().ToStdString(),tab_val_recon);
  //tableau de config  manip
  modif_tab_val("NXMAX",editNXMAX->GetValue().ToStdString(),tab_val_manip);

  modif_tab_val("VXMIN",editVxmin->GetValue().ToStdString(),tab_val_manip);
  modif_tab_val("VYMIN",editVymin->GetValue().ToStdString(),tab_val_manip);
  modif_tab_val("VXMAX",editVxmax->GetValue().ToStdString(),tab_val_manip);
  modif_tab_val("VYMAX",editVymax->GetValue().ToStdString(),tab_val_manip);

  modif_tab_val("CIRCLE_CX",editCX->GetValue().ToStdString(),tab_val_manip);
  modif_tab_val("CIRCLE_CY",editCY->GetValue().ToStdString(),tab_val_manip);


  sav_val(chemin_recon,tab_val_recon);
  sav_val(chemin_config_manip,tab_val_manip);
  sav_val(chemin_config_GUI,tab_val_gui_conf);
}

//extraire la valeur d un tokan dans le fichier situé a chemin_fic
/*float Gui_TomoFrame::extract_val(std::string token,  std::string chemin_fic)
{
    ifstream fichier(chemin_fic.c_str(), ios::in);  // on ouvre en lecture
    string ligne,motcle,valeurMot,separ=" ";
    float valeur=0;
    vector<std::string> tokens;

    if(fichier)  // si l'ouverture a fonctionné
    {
        while(!fichier.eof()){
            getline(fichier,ligne,'\n');//extrait chaque ligne du fichier (séparateur=retour chariot)
            if(ligne[0]!='#')//si pas ligne de commentaire
            tokens.push_back(ligne);
        }
    }
    else
        cerr << "Impossible d'ouvrir le fichier !"<< chemin_fic<< endl;

    int nb_tok=tokens.size();
    for(int cpt=0;cpt<nb_tok;cpt++){
        ligne=tokens[cpt];
        if(ligne!=""){
            int pos_separ=ligne.find(separ);
            int long_separ=separ.length();
            motcle = ligne.substr(0, pos_separ);//sbstr(pos_debut,pos_fin)
            if(motcle==token){
            valeurMot=ligne.substr(pos_separ+long_separ,ligne.size()-(motcle.size()+long_separ));
            cout<<motcle<<"="<<valeurMot<<endl;
            valeur=atof(valeurMot.c_str());
            }
        }
    }
    fichier.close();
    return valeur;
}*/

string Gui_TomoFrame::extract_string(std::string token,  std::string chemin_fic)
{
    ifstream fichier(chemin_fic.c_str(), ios::in);  // on ouvre en lecture
    string ligne,motcle,valeurMot,separ=" ";
    string valeur;
    vector<std::string> tokens;

    if(fichier)  // si l'ouverture a fonctionné
    {
        while(!fichier.eof()){
            getline(fichier,ligne,'\n');//extrait chaque ligne du fichier (séparateur=retour chariot)
            if(ligne[0]!='#')//si pas ligne de commentaire
            tokens.push_back(ligne);
        }
    }
    else
        cerr << "Impossible d'ouvrir le fichier !"<< chemin_fic<< endl;

    int nb_tok=tokens.size();
    for(int cpt=0;cpt<nb_tok;cpt++){
        ligne=tokens[cpt];
        if(ligne!=""){
            int pos_separ=ligne.find(separ);
            int long_separ=separ.length();
            motcle = ligne.substr(0, pos_separ);//sbstr(pos_debut,pos_fin)
            if(motcle==token){
            valeurMot=ligne.substr(pos_separ+long_separ,ligne.size()-(motcle.size()+long_separ));
            cout<<motcle<<"="<<valeurMot<<endl;
            valeur=valeurMot.c_str();
            }
        }
    }
    if(valeur.empty())
        cout<<"mot_clé "<<token<<" inexistant dans le fichier "<<chemin_fic<<endl;
    fichier.close();
    return valeur;
}

///charger en mémoire le tableau de mot-clés à partir d'un fichier de config.
void Gui_TomoFrame::init_tab_val(string chemin_fic, vector<string> &tab_val)
{     //effacer l'ancienne version->il faudrait initialiser dans le corps du programme!
    while (!tab_val.empty())
        tab_val.pop_back();

     string ligne;
     ifstream flux_fichier(chemin_fic.c_str(), ios::in);//ouverture en lecture ecriture

     if(flux_fichier)  // si l'ouverture a réussi
     {    while(!flux_fichier.eof()){
            getline(flux_fichier,ligne,'\n');//extrait chaque ligne du fichier (séparateur=retour chariot)
            tab_val.push_back(ligne);
            }
        }
     else  // sinon erreur
            cerr << "Erreur d'ouverture du fichier"<<chemin_fic << endl;
    flux_fichier.close();
}


//modifier en mémoire la valeur des token
void Gui_TomoFrame::modif_tab_val(string token,string string_valeur_token,vector<string> &tab_val)
{
          std::vector<std::string> fichier_new;
          string ligne,motcle,string_val_token,nlle_ligne,separ=" ";

          for(size_t cpt=0;cpt<tab_val.size();cpt++){
            ligne=tab_val[cpt];

           if(ligne!='#')//si pas ligne de commentaire (premier caractère ligne!=#)
            {
                int pos_separ=ligne.find(separ);//trouver l'espace (séparateur mot-clé valeur)
               //cout<<"position de l'espace="<<pos_separ<<endl;
                int long_separ=separ.length();
                motcle = ligne.substr(0, pos_separ);//copie une portion de lachaine entre 0 et pos_separ=isole le mot-clé
                if(motcle==token){
                nlle_ligne=token+" "+string_valeur_token;
                //cout<<"nouvelle ligne"<<nlle_ligne<<endl;
                fichier_new.push_back(nlle_ligne);
                }
                else
                    fichier_new.push_back(ligne);//le motclé n'est pas le bon, on ecrit juste la ligne
            }
            else{///caratère=# :  ligne de commentaire
           fichier_new.push_back(ligne);
            }
          }
            //effacer l'ancienne version
            while (!tab_val.empty())
            {
                tab_val.pop_back();
            }
           //insérer la nouvelle version
            for(int cpt=0;cpt<fichier_new.size();cpt++){
                tab_val.push_back(fichier_new[cpt]);
           }
             //effacer le tampon
              for(int cpt=0;cpt<fichier_new.size();cpt++)
            {
                fichier_new.pop_back();
            }
}

void Gui_TomoFrame::sav_val(string chemin_fic,vector<string> &tab_val)
{
     ///on sauvegarde le fichier, au cas où...
     ifstream src(chemin_fic.c_str() ,ios::binary);
     ofstream dst(chemin_fic+"_SAV" ,ios::binary);
     dst<<src.rdbuf();
     src.close();
     dst.close();

    ofstream flux_ecriture(chemin_fic.c_str(), ios::out);  // ouverture en écriture
    for(int cpt=0;cpt<tab_val.size();cpt++){
    flux_ecriture<<tab_val[cpt]<<endl;//ecriture dans le fichier
    }
    flux_ecriture.close();  //  fermer

}




