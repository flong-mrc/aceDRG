
/*
 * File:   periodcTable.cpp
 * Author: flong
 *
 * Created on July 29, 2012, 11:22 PM
 */


#include "periodicTable.h"

namespace LIBMOL
{
    PeriodicTable::PeriodicTable()
    {
        matTypes[1]     = "H";                      // H
        matTypes[2]     = "Non-Metal";              // C, N, O, P, S, Se
        matTypes[3]     = "Alkali-Metals";          // Li, Na, K, Rb, Cs, Fr
        matTypes[4]     = "Alkaline-earth-Metals";  // Be, Mg, Ca, Sr, Ba, Ra
        matTypes[5]     = "Transition-Metals";      // Fe, Cu, Au, Mn, Ni etc
        matTypes[6]     = "Other-Metal";            // Al, Ga, In, Sn, Tl, Pb, Bi, etc.
        matTypes[7]     = "Semimetallics";          // B, Si, Ge, As, Sb, Te, Po
        matTypes[8]     = "Halogens";               // F, Cl, Br, I, At
        matTypes[9]     = "Rare-Earth";             // La, Ce etc.
        matTypes[10]    = "Inert-Gases";            // He, Ne, Ar, Kr, Xe, Rn

        /* For each element, the following properties are listed:
         *  1. row #   in the periodic table
         *  2. group # (general groups include the transition metal
         *  3. material type, such as metal, non-metal, rare-earth etc
         *  4. atomic number., (indication of the weight)
         *  5. default valence.
         *  6. default oxidation number. If it is 0, then it has no default value and
         *     needs to to be decided on-fly
         *
         *  Unknown  = -1
         */

        // Some of covalent radii for non-organic elements are from "Dalton Trans., 2008, 2832–2838,
        // Covalent radii revisited, by Beatriz Cordero at el"
        //
        // Hydrogen
        elements["H"]["row"]     = 1;
        elements["H"]["group"]   = 1;
        elements["H"]["matType"] = 1;
        elements["H"]["atomNum"] = 1;
        elements["H"]["val"]     = 1;
        elemProps["H"]["vdw"]    = 1.20;
        elemProps["H"]["cova"]   = 0.32;
        elemProps["H"]["ionM-"]   = 0.23;
        elemProps["H"]["ionM+"]   = 0.23;

        elements["D"]["row"]     = 1;
        elements["D"]["group"]   = 1;
        elements["D"]["matType"] = 1;
        elements["D"]["atomNum"] = 1;
        elements["D"]["val"]     = 1;
        elemProps["D"]["vdw"]    = 1.20;
        elemProps["D"]["cova"]   = 0.32;
        elemProps["D"]["ionM-"]   = 0.23;
        elemProps["D"]["ionM+"]   = 0.23;

        // Organic set or Non-metal
        elements["C"]["row"]     = 2;
        elements["C"]["group"]   = 14;
        elements["C"]["matType"] = 2;
        elements["C"]["atomNum"] = 6;
        elements["C"]["val"]     = 4;
        elemProps["C"]["vdw"]    = 1.70;
        elemProps["C"]["cova"]   = 0.85;    // 0.76;
        elemProps["C"]["ionM-"]   = 0.30;


        elements["N"]["row"]       = 2;
        elements["N"]["group"]     = 15;
        elements["N"]["matType"]   = 2;
        elements["N"]["atomNum"]   = 7;
        elements["N"]["val"]       = 3;
        elemProps["N"]["vdw"]      = 1.55;
        elemProps["N"]["cova"]     = 0.71;
        elemProps["N"]["ionM-"]    = 1.32;
        elemProps["N"]["ionM+"]    = 0.30;



        elements["O"]["row"]       = 2;
        elements["O"]["group"]     = 16;
        elements["O"]["matType"]   = 2;
        elements["O"]["atomNum"]   = 8;
        elements["O"]["val"]       = 2;
        elemProps["O"]["vdw"]      = 1.52;
        elemProps["O"]["cova"]     = 0.68;
        elemProps["O"]["ionM-"]    = 1.26;

        elements["P"]["row"]       = 3;
        elements["P"]["group"]     = 15;
        elements["P"]["matType"]   = 2;
        elements["P"]["atomNum"]   = 15;
        elements["P"]["val"]       = 5;
        elemProps["P"]["vdw"]      = 1.52;
        elemProps["P"]["vdw"]      = 1.80;
        elemProps["P"]["cova"]     = 1.05;
        elemProps["P"]["ionM+"]   = 0.58;
        //elemProps["P5+"]["cova"]   = 0.52;

        elements["S"]["row"]       = 3;
        elements["S"]["group"]     = 16;
        elements["S"]["matType"]   = 2;
        elements["S"]["atomNum"]   = 16;
        elements["S"]["val"]       = 6;
        elemProps["S"]["vdw"]      = 1.80;
        elemProps["S"]["vdw"]      = 1.80;
        elemProps["S"]["cova"]     = 1.02;
        elemProps["S"]["ionM-"]    = 1.70;
        elemProps["S"]["ionM+"]    = 0.51;
        //elemProps["S6+"]["cova"]   = 0.43;

        elements["Se"]["row"]      = 4;
        elements["Se"]["group"]    = 16;
        elements["Se"]["matType"]  = 2;
        elements["Se"]["atomNum"]  = 34;
        elements["Se"]["val"]      = 6;
        elemProps["Se"]["vdw"]     = 1.90;
        elemProps["Se"]["cova"]    = 1.16;
        elemProps["Se"]["ionM-"]   = 1.84;
        elemProps["Se"]["ionM+"]   = 0.64;

        //elemProps["Se6+"]["cova"]   = 0.56;
        elements["B"]["row"]        = 2;
        elements["B"]["group"]      = 13;
        elements["B"]["matType"]    = 7;
        elements["B"]["atomNum"]    = 5;
        elements["B"]["val"]        = 3;
        //elements["B"]["vdw"]        = 1.85;
        //elements["B"]["cova"]       = 0.83;
        elemProps["B"]["vdw"]       = 1.85;
        elemProps["B"]["cova"]      = 0.83;
        elemProps["B"]["ionM+"]     = 0.41;


        // Halogens
        elements["F"]["row"]     = 2;
        elements["F"]["group"]   = 17;
        elements["F"]["matType"] = 8;
        elements["F"]["atomNum"] = 9;
        elements["F"]["val"]     = 1;
        //elements["F"]["vdw"]     = 1.47;
        //elements["F"]["cova"]    = 0.64;
        elemProps["F"]["vdw"]    = 1.47;
        elemProps["F"]["cova"]   = 0.57;
        elemProps["F"]["ionM-"] = 1.19;
        elemProps["F"]["ionM+"] = 0.22;

        elements["Cl"]["row"]     = 3;
        elements["Cl"]["group"]   = 17;
        elements["Cl"]["matType"] = 8;
        elements["Cl"]["atomNum"] = 17;
        elements["Cl"]["val"]     = 1;
        //elements["Cl"]["vdw"]     = 1.75;
        //elements["Cl"]["cova"]    = 0.99;
        elemProps["Cl"]["vdw"]    = 1.75;
        elemProps["Cl"]["cova"]    = 0.99;
        elemProps["Cl"]["ionM-"]  = 1.67;

        elements["Br"]["row"]     = 4;
        elements["Br"]["group"]   = 17;
        elements["Br"]["matType"] = 8;
        elements["Br"]["atomNum"] = 35;
        elements["Br"]["val"]     = 1;
        //elements["Br"]["vdw"]     = 1.85;
        //elements["Br"]["cova"]    = 1.21;
        elemProps["Br"]["vdw"]    = 1.85;
        elemProps["Br"]["cova"]    = 1.14;
        elemProps["Br"]["ionM-"]  = 1.96;

        elements["I"]["row"]      = 5;
        elements["I"]["group"]    = 17;
        elements["I"]["matType"]  = 8;
        elements["I"]["atomNum"]  = 53;
        elements["I"]["val"]      = 1;
        //elements["I"]["vdw"]      = 1.98;
        //elements["I"]["cova"]     = 1.40;
        elemProps["I"]["vdw"]     = 1.98;
        elemProps["I"]["cova"]     = 1.33;
        elemProps["I"]["ionM-"]   = 2.06;

        elements["At"]["row"]      = 6;
        elements["At"]["group"]    = 17;
        elements["At"]["matType"]  = 8;
        elements["At"]["atomNum"]  = 85;
        elements["At"]["val"]      = 1;
        //elements["At"]["vdw"]      = 2.00;
        //elements["At"]["cova"]     = 1.50;
        elemProps["At"]["vdw"]     = 2.00;
        elemProps["At"]["cova"]    = 1.50;
        elemProps["At"]["ionM-"]  = 0.76;

        // Alkali-Metals
        elements["Li"]["row"]      = 2;
        elements["Li"]["group"]    = 1;
        elements["Li"]["matType"]  = 3;
        elements["Li"]["atomNum"]  = 3;
        elements["Li"]["val"]      = 1;
        //elements["Li"]["vdw"]      = 1.82;
        //elements["Li"]["cova"]     = 0.68;
        elemProps["Li"]["vdw"]     = 1.82;
        elemProps["Li"]["cova"]    = 1.28;
        elemProps["Li"]["ionM+"]  = 0.73;

        elements["Na"]["row"]      = 3;
        elements["Na"]["group"]    = 1;
        elements["Na"]["matType"]  = 3;
        elements["Na"]["atomNum"]  = 11;
        elements["Na"]["val"]      = 1;
        //elements["Na"]["vdw"]      = 2.27;
        //elements["Na"]["cova"]     = 0.97;
        elemProps["Na"]["vdw"]     = 2.27;
        elemProps["Na"]["cova"]    = 1.13;
        elemProps["Na"]["ionM+"]   = 1.13;

        elements["K"]["row"]       = 4;
        elements["K"]["group"]     = 1;
        elements["K"]["matType"]   = 3;
        elements["K"]["atomNum"]   = 19;
        elements["K"]["val"]       = 1;
        //elements["K"]["vdw"]       = 2.75;
        //elements["K"]["cova"]      = 1.33;
        elemProps["K"]["vdw"]      = 2.75;
        elemProps["K"]["cova"]     = 1.23;
        elemProps["K"]["ionM+"]   = 1.52;

        elements["Rb"]["row"]      = 5;
        elements["Rb"]["group"]    = 1;
        elements["Rb"]["matType"]  = 3;
        elements["Rb"]["atomNum"]  = 37;
        elements["Rb"]["val"]      = 1;
        //elements["Rb"]["vdw"]      = 2.00;
        //elements["Rb"]["cova"]     = 1.47;
        elemProps["Rb"]["vdw"]     = 2.00;
        elemProps["Rb"]["cova"]    = 1.48;
        elemProps["Rb"]["ionM+"]  = 1.48;

        elements["Cs"]["row"]      = 6;
        elements["Cs"]["group"]    = 1;
        elements["Cs"]["matType"]  = 3;
        elements["Cs"]["atomNum"]  = 55;
        elements["Cs"]["val"]      = 1;
        elements["Cs"]["vdw"]      = 2.00;
        //elements["Cs"]["cova"]     = 1.67;
        //elemProps["Cs"]["vdw"]     = 2.00;
        elemProps["Cs"]["cova"]    = 1.67;
        elemProps["Cs"]["ionM+"]    = 1.81;

        elements["Fr"]["row"]      = 7;
        elements["Fr"]["group"]    = 1;
        elements["Fr"]["matType"]  = 3;
        elements["Fr"]["atomNum"]  = 87;
        elements["Fr"]["val"]      = 1;
        //elements["Fr"]["vdw"]      = 2.00;
        //elements["Fr"]["cova"]     = 1.50;
        elemProps["Fr"]["vdw"]     = 2.00;    // need to find out it
        elemProps["Fr"]["cova"]    = 1.50;
        elemProps["Fr"]["ionM+"]    = 1.94;

        // Alkaline-earth-Metals
        elements["Be"]["row"]      = 2;
        elements["Be"]["group"]    = 2;
        elements["Be"]["matType"]  = 4;
        elements["Be"]["atomNum"]  = 4;
        elements["Be"]["val"]      = 2;
        //elements["Be"]["vdw"]      = 2.00;
        //elements["Be"]["cova"]     = 0.35;
        elemProps["Be"]["vdw"]     = 2.00;
        elemProps["Be"]["cova"]    = 0.96;
        elemProps["Be"]["ionM+"]  = 0.41;

        elements["Mg"]["row"]      = 3;
        elements["Mg"]["group"]    = 2;
        elements["Mg"]["matType"]  = 4;
        elements["Mg"]["atomNum"]  = 12;
        elements["Mg"]["val"]      = 2;
        //elements["Mg"]["vdw"]      = 1.73;
        //elements["Mg"]["cova"]     = 1.10;
        elemProps["Mg"]["vdw"]     = 1.73;
        elemProps["Mg"]["cova"]     = 1.41;
        elemProps["Mg"]["ionM+"]   = 0.71;

        elements["Ca"]["row"]      = 4;
        elements["Ca"]["group"]    = 2;
        elements["Ca"]["matType"]  = 4;
        elements["Ca"]["atomNum"]  = 20;
        elements["Ca"]["val"]      = 2;
        //elements["Ca"]["vdw"]      = 2.00;
        //elements["Ca"]["cova"]     = 0.99;
        elemProps["Ca"]["vdw"]     = 2.00;
        elemProps["Ca"]["cova"]    = 1.74;
        elemProps["Ca"]["ionM+"]  = 1.14;

        elements["Sr"]["row"]      = 5;
        elements["Sr"]["group"]    = 2;
        elements["Sr"]["matType"]  = 4;
        elements["Sr"]["atomNum"]  = 38;
        elements["Sr"]["val"]      = 2;
        //elements["Sr"]["vdw"]      = 2.19;
        //elements["Sr"]["cova"]     = 1.12;
        elemProps["Sr"]["vdw"]     = 2.19;
        // elemProps["Sr"]["cova"]    = 1.92;
        elemProps["Sr"]["cova"]    = 1.85;
        elemProps["Sr"]["ionM+"]    = 1.32;

        elements["Ba"]["row"]      = 6;
        elements["Ba"]["group"]    = 2;
        elements["Ba"]["matType"]  = 4;
        elements["Ba"]["atomNum"]  = 56;
        elements["Ba"]["val"]      = 2;
        //elements["Ba"]["vdw"]      = 2.53;
        //elements["Ba"]["cova"]     = 1.34;
        elemProps["Ba"]["vdw"]     = 2.53;
        elemProps["Ba"]["cova"]    = 1.98;
        elemProps["Ba"]["ionM+"]  = 1.49;

        elements["Ra"]["row"]      = 7;
        elements["Ra"]["group"]    = 2;
        elements["Ra"]["matType"]  = 4;
        elements["Ra"]["atomNum"]  = 88;
        elements["Ra"]["val"]      = 2;
        //elements["Ra"]["vdw"]      = 2.15;
        //elements["Ra"]["cova"]     = 1.90;
        elemProps["Ra"]["vdw"]     = 2.15;
        elemProps["Ra"]["cova"]    = 1.90;
        elemProps["Ra"]["ionM+"]    = 1.62;


        // Transition-Metals

        elements["Sc"]["row"]      = 4;
        elements["Sc"]["group"]    = 3;
        elements["Sc"]["matType"]  = 5;
        elements["Sc"]["atomNum"]  = 21;
        elements["Sc"]["val"]      = 3;
        //elements["Sc"]["vdw"]      = 1.60;
        //elements["Sc"]["cova"]     = 1.44;
        elemProps["Sc"]["vdw"]     = 1.60;
        elemProps["Sc"]["cova"]    = 1.70;
        elemProps["Sc"]["ionM+"]  = 0.885;

        elements["Y"]["row"]       = 5;
        elements["Y"]["group"]     = 3;
        elements["Y"]["matType"]   = 5;
        elements["Y"]["atomNum"]   = 39;
        elements["Y"]["val"]       = 3;
        //elements["Y"]["vdw"]       = 1.80;
        //elements["Y"]["cova"]      = 1.78;
        elemProps["Y"]["vdw"]      = 1.80;
        elemProps["Y"]["cova"]     = 1.90;
        elemProps["Y"]["ionM+"]   = 1.04;

        elements["Ti"]["row"]      = 4;
        elements["Ti"]["group"]    = 4;
        elements["Ti"]["matType"]  = 5;
        elements["Ti"]["atomNum"]  = 22;
        elements["Ti"]["val"]      = 4;
        //elements["Ti"]["vdw"]      = 1.40;
        //elements["Ti"]["cova"]     = 1.47;
        elemProps["Ti"]["vdw"]     = 1.40;
        elemProps["Ti"]["cova"]    = 1.60;
        elemProps["Ti"]["ionM+"]  = 0.86;
        //elemProps["Ti3+"]["cova"]  = 0.67;
        //elemProps["Ti4+"]["cova"]  = 0.61;

        elements["Zr"]["row"]      = 5;
        elements["Zr"]["group"]    = 4;
        elements["Zr"]["matType"]  = 5;
        elements["Zr"]["atomNum"]  = 40;
        elements["Zr"]["val"]      = 4;
        //elements["Zr"]["vdw"]      = 1.55;
        //elements["Zr"]["cova"]     = 1.56;
        elemProps["Zr"]["vdw"]     = 1.55;
        elemProps["Zr"]["cova"]    = 1.75;
        elemProps["Zr"]["ionM+"]  = 0.72;

        elements["Hf"]["row"]      = 6;
        elements["Hf"]["group"]    = 4;
        elements["Hf"]["matType"]  = 5;
        elements["Hf"]["atomNum"]  = 72;
        elements["Hf"]["val"]      = 4;
        //elements["Hf"]["vdw"]      = 1.55;
        //elements["Hf"]["cova"]     = 1.57;
        elemProps["Hf"]["vdw"]     = 1.55;
        elemProps["Hf"]["cova"]     = 1.75;
        elemProps["Hf"]["ionM+"]   = 0.85;


        elements["Rf"]["row"]      = 7;
        elements["Rf"]["group"]    = 4;
        elements["Rf"]["matType"]  = 5;
        elements["Rf"]["atomNum"]  = 104;
        elements["Rf"]["val"]      = 4;
        elemProps["Rf"]["vdw"]     = 1.60;    // need to find it out
        elemProps["Rf"]["cova"]     = 1.21;
        elemProps["Rf"]["ionM+"]   = 0.85;

        elements["V"]["row"]       = 4;
        elements["V"]["group"]     = 5;
        elements["V"]["matType"]   = 5;
        elements["V"]["atomNum"]   = 23;
        elements["V"]["val"]      = 5;
        //elements["V"]["vdw"]      = 1.35;
        //elements["V"]["cova"]     = 1.33;
        elemProps["V"]["vdw"]     = 1.35;
        elemProps["V"]["cova"]    = 1.43;  //1.53;
        elemProps["V"]["ionM+"]   = 0.60;  // 0.68;
        //elemProps["V3+"]["cova"]    = 0.78;
        //elemProps["V4+"]["cova"]    = 0.72;

        elements["Nb"]["row"]      = 5;
        elements["Nb"]["group"]    = 5;
        elements["Nb"]["matType"]  = 5;
        elements["Nb"]["atomNum"]  = 41;
        elements["Nb"]["val"]      = 5;
        //elements["Nb"]["vdw"]      = 1.45;
        //elements["Nb"]["cova"]     = 1.48;
        elemProps["Nb"]["vdw"]     = 1.45;
        elemProps["Nb"]["cova"]    = 1.64;
        elemProps["Nb"]["ionM+"]  = 0.62;
        //elemProps["Nb4+"]["cova"]  = 0.82;
        //elemProps["Nb5+"]["cova"]  = 0.78;

        elements["Ta"]["row"]       = 6;
        elements["Ta"]["group"]     = 5;
        elements["Ta"]["matType"]   = 5;
        elements["Ta"]["atomNum"]   = 73;
        elements["Ta"]["val"]       = 5;
        //elements["Ta"]["vdw"]       = 1.45;
        //elements["Ta"]["cova"]      = 1.43;
        elemProps["Ta"]["vdw"]      = 1.45;
        elemProps["Ta"]["cova"]     = 1.70;
        elemProps["Ta"]["ionM+"]   = 0.78;
        //elemProps["Ta4+"]["cova"]   = 0.82;
        //elemProps["Ta5+"]["cova"]   = 0.78;


        elements["Db"]["row"]       = 7;
        elements["Db"]["group"]     = 5;
        elements["Db"]["matType"]   = 5;
        elements["Db"]["atomNum"]   = 105;
        elements["Db"]["val"]       = 5;
        //elements["Db"]["vdw"]       = 1.45;
        //elements["Db"]["cova"]      = 1.50;
        elemProps["Db"]["vdw"]      = 1.45;    // need to find it out
        elemProps["Db"]["cova"]     = 1.50;
        elemProps["Db"]["ionM+"]   = 0.78;


        elements["Cr"]["row"]       = 4;
        elements["Cr"]["group"]     = 6;
        elements["Cr"]["matType"]   = 5;
        elements["Cr"]["atomNum"]   = 24;
        elements["Cr"]["val"]       = 6;
        //elements["Cr"]["vdw"]       = 1.40;
        //elements["Cr"]["cova"]      = 1.35;
        elemProps["Cr"]["vdw"]      = 1.40;
        elemProps["Cr"]["cova"]     = 1.39;
        elemProps["Cr"]["ionM+"]    = 0.87;
        //elemProps["Cr3+"]["cova"]   = 0.76;
        //elemProps["Cr4+"]["cova"]   = 0.69;
        //elemProps["Cr5+"]["cova"]   = 0.63;
        //elemProps["Cr6+"]["cova"]   = 0.58;

        elements["Mo"]["row"]       = 5;
        elements["Mo"]["group"]     = 6;
        elements["Mo"]["matType"]   = 5;
        elements["Mo"]["atomNum"]   = 42;
        elements["Mo"]["val"]       = 6;
        //elements["Mo"]["vdw"]       = 1.45;
        //elements["Mo"]["cova"]      = 1.47;
        elemProps["Mo"]["vdw"]      = 1.45;
        elemProps["Mo"]["cova"]     = 1.29;     //1.39;
        elemProps["Mo"]["ionM+"]    = 0.69;
        //elemProps["Mo4+"]["cova"]   = 0.79;
        //elemProps["Mo5+"]["cova"]   = 0.75;
        //elemProps["Mo6+"]["cova"]   = 0.73;


        elements["W"]["row"]        = 6;
        elements["W"]["group"]      = 6;
        elements["W"]["matType"]    = 5;
        elements["W"]["atomNum"]    = 74;
        elements["W"]["val"]        = 6;
        //elements["W"]["vdw"]        = 1.85;
        //elements["W"]["cova"]       = 1.21;
        elemProps["W"]["vdw"]       = 1.35;
        elemProps["W"]["cova"]      = 1.62;
        elemProps["W"]["ionM+"]     = 0.80;
        //elemProps["W5+"]["cova"]   = 0.76;
        //elemProps["W6+"]["cova"]   = 0.74;

        elements["Sg"]["row"]       = 7;
        elements["Sg"]["group"]     = 6;
        elements["Sg"]["matType"]   = 5;
        elements["Sg"]["atomNum"]   = 106;
        elements["Sg"]["val"]       = 6;
        //elements["Sg"]["vdw"]       = 1.35;
        //elements["Sg"]["cova"]      = 1.50;
        elemProps["Sg"]["vdw"]      = 1.35;
        elemProps["Sg"]["cova"]     = 1.50;
        elemProps["Sg"]["ionM+"]    = 0.79;

        elements["Mn"]["row"]       = 4;
        elements["Mn"]["group"]     = 7;
        elements["Mn"]["matType"]   = 5;
        elements["Mn"]["atomNum"]   = 25;
        elements["Mn"]["val"]       = 7;
        //elements["Mn"]["vdw"]       = 1.40;
        //elements["Mn"]["cova"]      = 1.35;
        elemProps["Mn"]["vdw"]      = 1.40;
        elemProps["Mn"]["cova"]     = 1.39;
        elemProps["Mn"]["ionM+"]   = 0.81;
        //elemProps["Mn4+"]["cova"]   = 0.72;
        //elemProps["Mn5+"]["cova"]   = 0.67;
        //elemProps["Mn6+"]["cova"]   = 0.47;
        //elemProps["Mn7+"]["cova"]   = 0.40;

        elements["Tc"]["row"]       = 5;
        elements["Tc"]["group"]     = 7;
        elements["Tc"]["matType"]   = 5;
        elements["Tc"]["atomNum"]   = 43;
        elements["Tc"]["val"]       = 7;
        //elements["Tc"]["vdw"]       = 1.35;
        //elements["Tc"]["cova"]      = 1.35;
        elemProps["Tc"]["vdw"]      = 1.35;
        elemProps["Tc"]["cova"]     = 1.47;
        elemProps["Tc"]["ionM+"]    = 0.79;
        //elemProps["Tc5+"]["cova"]   = 0.64;
        //elemProps["Tc7+"]["cova"]   = 0.70;

        elements["Re"]["row"]       = 6;
        elements["Re"]["group"]     = 7;
        elements["Re"]["matType"]   = 5;
        elements["Re"]["atomNum"]   = 75;
        elements["Re"]["val"]       = 7;
        //elements["Re"]["vdw"]       = 1.35;
        //elements["Re"]["cova"]      = 1.35;
        elemProps["Re"]["vdw"]      = 1.35;
        elemProps["Re"]["cova"]     = 1.51;
        elemProps["Re"]["ionM+"]    = 0.77;
        //elemProps["Re5+"]["cova"]   = 0.72;
        //elemProps["Re6+"]["cova"]   = 0.69;
        //elemProps["Re7+"]["cova"]   = 0.67;

        elements["Bh"]["row"]       = 7;
        elements["Bh"]["group"]     = 7;
        elements["Bh"]["matType"]   = 5;
        elements["Bh"]["atomNum"]   = 107;
        elements["Bh"]["val"]       = 7;
        //elements["Bh"]["vdw"]       = 1.35;
        //elements["Bh"]["cova"]      = 1.50;
        elemProps["Bh"]["vdw"]      = 1.35;   // Need to find it out
        elemProps["Bh"]["cova"]     = 1.50;
        elemProps["Bh"]["ionM+"]    = 0.77;


        elements["Fe"]["row"]       = 4;
        elements["Fe"]["group"]     = 8;
        elements["Fe"]["matType"]   = 5;
        elements["Fe"]["atomNum"]   = 26;
        elements["Fe"]["val"]       = 6;
        //elements["Fe"]["vdw"]       = 1.40;
        //elements["Fe"]["cova"]      = 1.34;
        elemProps["Fe"]["vdw"]      = 1.40;
        elemProps["Fe"]["cova"]     = 1.32;
        elemProps["Fe"]["ionM+"]    = 0.75;
        //elemProps["Fe3+"]["cova"]   = 0.69;
        //elemProps["Fe4+"]["cova"]   = 0.725;
        //elemProps["Fe6+"]["cova"]   = 0.39;



        elements["Ru"]["row"]       = 5;
        elements["Ru"]["group"]     = 8;
        elements["Ru"]["matType"]   = 5;
        elements["Ru"]["atomNum"]   = 44;
        elements["Ru"]["val"]       = 8;
        //elements["Ru"]["vdw"]       = 1.30;
        //elements["Ru"]["cova"]      = 1.40;
        elemProps["Ru"]["vdw"]      = 1.30;
        elemProps["Ru"]["cova"]     = 1.46;
        elemProps["Ru"]["ionM+"]   = 0.82;
        //elemProps["Ru4+"]["cova"]   = 0.76;
        //elemProps["Ru5+"]["cova"]   = 0.705;
        //elemProps["Ru7+"]["cova"]   = 0.52;
        //elemProps["Ru8+"]["cova"]   = 0.50;

        elements["Os"]["row"]       = 6;
        elements["Os"]["group"]     = 8;
        elements["Os"]["matType"]   = 5;
        elements["Os"]["atomNum"]   = 76;
        elements["Os"]["val"]       = 8;
        //elements["Os"]["vdw"]       = 1.30;
        //elements["Os"]["cova"]      = 1.37;
        elemProps["Os"]["vdw"]      = 1.30;
        elemProps["Os"]["cova"]      = 1.44;
        elemProps["Os"]["ionM+"]   = 0.77;
        //elemProps["Os5+"]["cova"]   = 0.715;
        //elemProps["Os6+"]["cova"]   = 0.685;
        //elemProps["Os7+"]["cova"]   = 0.665;
        //elemProps["Os8+"]["cova"]   = 0.53;

        elements["Hs"]["row"]       = 7;
        elements["Hs"]["group"]     = 8;
        elements["Hs"]["matType"]   = 5;
        elements["Hs"]["atomNum"]   = 108;
        elements["Hs"]["val"]       = 8;
        //elements["Hs"]["vdw"]       = 1.30;
        //elements["Hs"]["cova"]      = 1.50;
        elemProps["Hs"]["vdw"]      = 1.30;   // Need to find it out
        elemProps["Hs"]["cova"]     = 1.50;
        elemProps["Hs"]["ionM+"]   = 0.805;

        elements["Co"]["row"]       = 4;
        elements["Co"]["group"]     = 9;
        elements["Co"]["matType"]   = 5;
        elements["Co"]["atomNum"]   = 27;
        elements["Co"]["val"]       = 5;
        //elements["Co"]["vdw"]       = 1.35;
        //elements["Co"]["cova"]      = 1.33;
        elemProps["Co"]["vdw"]      = 1.35;
        elemProps["Co"]["cova"]     = 1.33; // 1.26
        elemProps["Co"]["ionM+"]    = 0.885;
        //elemProps["Co3+"]["cova"]   = 0.75;
        //elemProps["Co4+"]["cova"]   = 0.67;


        elements["Rh"]["row"]       = 5;
        elements["Rh"]["group"]     = 9;
        elements["Rh"]["matType"]   = 5;
        elements["Rh"]["atomNum"]   = 45;
        //elements["Rh"]["val"]       = 6;
        //elements["Rh"]["vdw"]       = 1.35;
        elemProps["Rh"]["vdw"]      = 1.35;
        elemProps["Rh"]["cova"]     = 1.42;
        elemProps["Rh"]["ionM+"]   = 0.805;
        //elemProps["Rh4+"]["cova"]   = 0.74;
        //elemProps["Rh5+"]["cova"]   = 0.69;


        elements["Ir"]["row"]       = 6;
        elements["Ir"]["group"]     = 9;
        elements["Ir"]["matType"]   = 5;
        elements["Ir"]["atomNum"]   = 77;
        elements["Ir"]["val"]       = 8;
        //elements["Ir"]["vdw"]       = 1.35;
        //elements["Ir"]["cova"]      = 1.32;
        elemProps["Ir"]["vdw"]      = 1.35;
        elemProps["Ir"]["ionM+"]    = 0.82;
        elemProps["Ir"]["cova"]     = 1.41;
        //elemProps["Ir5+"]["cova"]   = 0.71;

        elements["Mt"]["row"]       = 7;
        elements["Mt"]["group"]     = 9;
        elements["Mt"]["matType"]   = 5;
        elements["Mt"]["atomNum"]   = 109;
        elements["Mt"]["val"]       = -1;
        //elements["Mt"]["vdw"]       = 1.35;
        //elements["Mt"]["cova"]      = 1.21;
        elemProps["Mt"]["vdw"]      = 1.35;   // Need to find it out
        elemProps["Mt"]["cova"]      = 1.21;
        elemProps["Mt"]["ionM+"]    = 0.82;

        elements["Ni"]["row"]       = 4;
        elements["Ni"]["group"]     = 10;
        elements["Ni"]["matType"]   = 5;
        elements["Ni"]["atomNum"]   = 28;
        elements["Ni"]["val"]       = 4;
        //elements["Ni"]["vdw"]       = 1.63;
        //elements["Ni"]["cova"]      = 1.50;
        elemProps["Ni"]["vdw"]      = 1.63;
        elemProps["Ni"]["cova"]     = 1.24;
        elemProps["Ni"]["ionM+"]    = 0.83;
        //elemProps["Ni3+"]["cova"]   = 0.70;
        //elemProps["Ni4+"]["cova"]   = 0.62;

        elements["Pd"]["row"]       = 5;
        elements["Pd"]["group"]     = 10;
        elements["Pd"]["matType"]   = 5;
        elements["Pd"]["atomNum"]   = 46;
        elements["Pd"]["val"]       = 4;
        //elements["Pd"]["vdw"]       = 1.63;
        //elements["Pd"]["cova"]      = 1.50;
        elemProps["Pd"]["vdw"]      = 1.63;
        elemProps["Pd"]["cova"]     =  1.50;  //1.39;
        elemProps["Pd"]["ionM+"]    = 0.73;
        //elemProps["Pd2+"]["cova"]    = 1.00;
        //elemProps["Pd3+"]["cova"]      = 0.90;
        //elemProps["Pd4+"]["cova"]      = 0.755;

        elements["Pt"]["row"]       = 6;
        elements["Pt"]["group"]     = 10;
        elements["Pt"]["matType"]   = 5;
        elements["Pt"]["atomNum"]   = 78;
        elements["Pt"]["val"]       = 6;
        //elements["Pt"]["vdw"]       = 1.72;
        //elements["Pt"]["cova"]      = 1.50;
        elemProps["Pt"]["vdw"]      = 1.72;
        elemProps["Pt"]["cova"]      = 1.60; // 1.36;
        elemProps["Pt"]["ionM+"]    = 0.94;
        //elemProps["Pt4+"]["cova"]    = 0.765;
        //elemProps["Pt5+"]["cova"]    = 0.71;

        elements["Ds"]["row"]       = 7;
        elements["Ds"]["group"]     = 10;
        elements["Ds"]["matType"]   = 5;
        elements["Ds"]["atomNum"]   = 110;
        elements["Ds"]["val"]       = -1;
        //elements["Ds"]["vdw"]       = 1.75;
        //elements["Ds"]["cova"]      = 1.50;
        elemProps["Ds"]["vdw"]      = 1.75;
        elemProps["Ds"]["cova"]      = 1.50;
        elemProps["Ds"]["ionM+"]    = 0.94;

        elements["Cu"]["row"]       = 4;
        elements["Cu"]["group"]     = 11;
        elements["Cu"]["matType"]   = 5;
        elements["Cu"]["atomNum"]   = 29;
        elements["Cu"]["val"]       = 2;
        //elements["Cu"]["vdw"]       = 1.40;
        //elements["Cu"]["cova"]      = 1.52;
        elemProps["Cu"]["vdw"]      = 1.40;
        elemProps["Cu"]["cova"]     = 1.32;
        elemProps["Cu"]["ionM+"]    = 0.91;
        //elemProps["Cu2+"]["cova"]    = 0.87;
        //elemProps["Cu3+"]["cova"]    = 0.68;

        elements["Ag"]["row"]       = 5;
        elements["Ag"]["group"]     = 11;
        elements["Ag"]["matType"]   = 5;
        elements["Ag"]["atomNum"]   = 47;
        elements["Ag"]["val"]       = 4;
        //elements["Ag"]["vdw"]       = 1.72;
        //elements["Ag"]["cova"]      = 1.59;
        elemProps["Ag"]["vdw"]      = 1.72;
        elemProps["Ag"]["cova"]     = 1.45;
        elemProps["Ag"]["ionM+"]    = 1.29;
        //elemProps["Ag2+"]["cova"]    = 1.08;
        //elemProps["Ag3+"]["cova"]    = 0.89;

        elements["Au"]["row"]       = 6;
        elements["Au"]["group"]     = 11;
        elements["Au"]["matType"]   = 5;
        elements["Au"]["atomNum"]   = 79;
        elements["Au"]["val"]       = 5;
        //elements["Au"]["vdw"]       = 1.66;
        //elements["Au"]["cova"]      = 1.50;
        elemProps["Au"]["vdw"]      = 1.66;
        elemProps["Au"]["cova"]     = 1.37;
        elemProps["Au"]["ionM+"]    = 1.51;
        //elemProps["Au3+"]["cova"]    = 0.99;
        //elemProps["Au5+"]["cova"]    = 0.71;

        elements["Rg"]["row"]       = 7;
        elements["Rg"]["group"]     = 11;
        elements["Rg"]["matType"]   = 5;
        elements["Rg"]["atomNum"]   = 111;
        elements["Rg"]["val"]       = -1;
        //elements["Rg"]["vdw"]       = 1.66;
        //elements["Rg"]["cova"]      = 1.50;
        elemProps["Rg"]["vdw"]      = 1.66;   // Need to find it out
        elemProps["Rg"]["cova"]     = 1.50;
        elemProps["Rg"]["ionM+"]    = 0.94;

        elements["Zn"]["row"]       = 4;
        elements["Zn"]["group"]     = 12;
        elements["Zn"]["matType"]   = 5;
        elements["Zn"]["atomNum"]   = 30;
        elements["Zn"]["val"]       = 2;
        //elements["Zn"]["vdw"]       = 1.39;
        //elements["Zn"]["cova"]      = 1.45;
        elemProps["Zn"]["vdw"]      = 1.39;
        elemProps["Zn"]["cova"]     = 1.18;
        elemProps["Zn"]["ionM+"]    = 0.88;

        elements["Cd"]["row"]       = 5;
        elements["Cd"]["group"]     = 12;
        elements["Cd"]["matType"]   = 5;
        elements["Cd"]["atomNum"]   = 48;
        elements["Cd"]["val"]       = 2;
        //elements["Cd"]["vdw"]       = 1.58;
        //elements["Cd"]["cova"]      = 1.69;
        elemProps["Cd"]["vdw"]      = 1.58;
        elemProps["Cd"]["cova"]     = 1.69; //     1.44;
        elemProps["Cd"]["ionM+"]    = 1.09;

        elements["Hg"]["row"]       = 6;
        elements["Hg"]["group"]     = 12;
        elements["Hg"]["matType"]   = 5;
        elements["Hg"]["atomNum"]   = 80;
        elements["Hg"]["val"]       = 2;
        //elements["Hg"]["vdw"]       = 155;
        //elements["Hg"]["cova"]      = 1.70;
        elemProps["Hg"]["vdw"]      = 1.55;
        elemProps["Hg"]["cova"]     = 1.70; //   1.32;
        elemProps["Hg"]["ionM+"]   = 1.33;
        //elemProps["Hg2+"]["cova"]   = 1.16;


        elements["Cn"]["row"]       = 7;
        elements["Cn"]["group"]     = 12;
        elements["Cn"]["matType"]   = 5;
        elements["Cn"]["atomNum"]   = 112;
        elements["Cn"]["val"]       = -1;
        elemProps["Cn"]["vdw"]      = 1.58;   // Need to find it out
        elemProps["Cn"]["cova"]     = 1.58;   // Need to find it out

        // Other-Metal
        elements["Al"]["row"]       = 3;
        elements["Al"]["group"]     = 13;
        elements["Al"]["matType"]   = 6;
        elements["Al"]["atomNum"]   = 13;
        elements["Al"]["val"]       = 3;
        //elements["Al"]["vdw"]       = 1.25;
        //elements["Al"]["cova"]      = 1.35;
        elemProps["Al"]["vdw"]      = 1.25;
        elemProps["Al"]["cova"]     = 1.21;
        elemProps["Al"]["ionM+"]    = 0.675;

        elements["Ga"]["row"]       = 4;
        elements["Ga"]["group"]     = 13;
        elements["Ga"]["matType"]   = 6;
        elements["Ga"]["atomNum"]   = 31;
        elements["Ga"]["val"]       = 3;
        //elements["Ga"]["vdw"]       = 1.87;
        //elements["Ga"]["cova"]      = 1.22;
        elemProps["Ga"]["vdw"]      = 1.87;
        elemProps["Ga"]["cova"]     = 1.22;
        elemProps["Ga"]["ionM+"]    = 0.76;

        elements["In"]["row"]       = 5;
        elements["In"]["group"]     = 13;
        elements["In"]["matType"]   = 6;
        elements["In"]["atomNum"]   = 49;
        elements["In"]["val"]       = 3;
        //elements["In"]["vdw"]       = 1.93;
        //elements["In"]["cova"]      = 1.63;
        elemProps["In"]["vdw"]      = 1.93;
        elemProps["In"]["cova"]     = 1.42;
        elemProps["In"]["ionM+"]    = 0.94;

        elements["Ti"]["row"]       = 6;
        elements["Ti"]["group"]     = 13;
        elements["Ti"]["matType"]   = 6;
        elements["Ti"]["atomNum"]   = 81;
        elements["Ti"]["val"]       = 3;
        //elements["Ti"]["vdw"]       = 1.40;
        //elements["Ti"]["cova"]      = 1.47;
        elemProps["Ti"]["vdw"]      = 1.40;
        elemProps["Ti"]["cova"]     = 1.60;
        elemProps["Ti"]["ionM+"]     = 1.00;
        //elemProps["Ti3+"]["cova"]     = 0.81;
        //elemProps["Ti4+"]["cova"]     = 0.745;

        elements["Sn"]["row"]       = 5;
        elements["Sn"]["group"]     = 14;
        elements["Sn"]["matType"]   = 6;
        elements["Sn"]["atomNum"]   = 50;
        elements["Sn"]["val"]       = 4;
        //elements["Sn"]["vdw"]       = 2.17;
        //elements["Sn"]["cova"]      = 1.46;
        elemProps["Sn"]["vdw"]      = 2.17;
        elemProps["Sn"]["cova"]     = 1.39;
        elemProps["Sn"]["ionM+"]     = 0.83;

        elements["Pb"]["row"]       = 6;
        elements["Pb"]["group"]     = 14;
        elements["Pb"]["matType"]   = 6;
        elements["Pb"]["atomNum"]   = 82;
        elements["Pb"]["val"]       = 4;
        //elements["Pb"]["vdw"]       = 2.02;
        //elements["Pb"]["cova"]      = 1.54;
        elemProps["Pb"]["vdw"]      = 2.02;
        elemProps["Pb"]["cova"]     = 1.46;
        elemProps["Pb"]["ionM+"]     = 1.33;
        //elemProps["Pb+4"]["cova"]     = 0.915;

        elements["Bi"]["row"]       = 6;
        elements["Bi"]["group"]     = 15;
        elements["Bi"]["matType"]   = 6;
        elements["Bi"]["atomNum"]   = 83;
        elements["Bi"]["val"]       = 5;
        //elements["Bi"]["vdw"]       = 1.80;
        //elements["Bi"]["cova"]      = 1.54;
        elemProps["Bi"]["vdw"]      = 1.80;
        elemProps["Bi"]["cova"]     = 1.42;  //1.48;
        elemProps["Bi"]["ionM+"]    = 1.07;
        // Semimetallics


        elements["Si"]["row"]       = 3;
        elements["Si"]["group"]     = 14;
        elements["Si"]["matType"]   = 7;
        elements["Si"]["atomNum"]   = 14;
        elements["Si"]["val"]       = 4;
        //elements["Si"]["vdw"]       = 2.10;
        //elements["Si"]["cova"]      = 1.20;
        elemProps["Si"]["vdw"]      = 2.10;
        elemProps["Si"]["cova"]     = 1.20;
        elemProps["Si"]["ionM+"]    = 0.54;

        elements["Ge"]["row"]       = 4;
        elements["Ge"]["group"]     = 14;
        elements["Ge"]["matType"]   = 7;
        elements["Ge"]["atomNum"]   = 32;
        elements["Ge"]["val"]       = 4;
        //elements["Ge"]["vdw"]       = 2.10;
        //elements["Ge"]["cova"]      = 1.17;
        elemProps["Ge"]["vdw"]      = 2.10;
        elemProps["Ge"]["cova"]     = 1.20;
        elemProps["Ge"]["ionM+"]   = 0.87;
        //elemProps["Ge4+"]["cova"]   = 0.67;

        elements["As"]["row"]       = 4;
        elements["As"]["group"]     = 15;
        elements["As"]["matType"]   = 7;
        elements["As"]["atomNum"]   = 33;
        elements["As"]["val"]       = 5;
        //elements["As"]["vdw"]       = 1.85;
        //elements["As"]["cova"]      = 1.21;
        elemProps["As"]["vdw"]      = 1.85;
        elemProps["As"]["cova"]     = 1.21;  // 1.11;
        elemProps["As"]["ionM+"]    = 0.72;
        //elemProps["As5+"]["cova"]    = 0.60;

        elements["Sb"]["row"]       = 5;
        elements["Sb"]["group"]     = 15;
        elements["Sb"]["matType"]   = 7;
        elements["Sb"]["atomNum"]   = 51;
        elements["Sb"]["val"]       = 5;
        //elements["Sb"]["vdw"]       = 1.85;
        //elements["Sb"]["cova"]      = 1.21;
        elemProps["Sb"]["vdw"]      = 1.80;
        elemProps["Sb"]["cova"]     = 1.39;
        elemProps["Sb"]["ionM+"]    = 0.90;
        //elemProps["Sb5+"]["cova"]    = 0.74;

        elements["Te"]["row"]       = 5;
        elements["Te"]["group"]     = 16;
        elements["Te"]["matType"]   = 7;
        elements["Te"]["atomNum"]   = 52;
        elements["Te"]["val"]       = 6;
        //elements["Te"]["vdw"]       = 2.06;
        //elements["Te"]["cova"]      = 1.47;
        elemProps["Te"]["vdw"]      = 2.06;
        elemProps["Te"]["cova"]     = 1.38;
        elemProps["Te"]["ionM-"]    = 2.07;
        elemProps["Te"]["ionM+"]    = 1.11;
        //elemProps["Te5+"]["cova"]    = 0.70;

        elements["Po"]["row"]       = 6;
        elements["Po"]["group"]     = 16;
        elements["Po"]["matType"]   = 7;
        elements["Po"]["atomNum"]   = 84;
        elements["Po"]["val"]       = 6;
        //elements["Po"]["vdw"]       = 1.85;
        //elements["Po"]["cova"]      = 1.60;
        elemProps["Po"]["vdw"]      = 2.00;
        elemProps["Po"]["cova"]     = 1.40;
        elemProps["Po"]["ionM+"]    = 1.08;
        //elemProps["Po6+"]["cova"]    = 0.81;

        // Rare earth metals

        elements["La"]["row"]       = 6;
        elements["La"]["group"]     = 3;
        elements["La"]["matType"]   = 9;
        elements["La"]["atomNum"]   = 57;
        elements["La"]["val"]       = 3;
        //elements["Po"]["vdw"]       = 1.85;
        //elements["Po"]["cova"]      = 1.60;
        elemProps["La"]["vdw"]      = 1.87;  //atom radius
        elemProps["La"]["cova"]     = 2.07;
        elemProps["La"]["ionM+"]    = 1.08;


        elements["Ce"]["row"]       = 6;
        elements["Ce"]["group"]     = 3;
        elements["Ce"]["matType"]   = 9;
        elements["Ce"]["atomNum"]   = 58;
        elements["Ce"]["val"]       = 4;
        //elements["Po"]["vdw"]       = 1.85;
        //elements["Po"]["cova"]      = 1.60;
        elemProps["Ce"]["vdw"]      = 1.82;  //atom radius
        elemProps["Ce"]["cova"]     = 2.04;
        elemProps["Ce"]["ionM+"]    = 1.08;


        elements["Pr"]["row"]       = 6;
        elements["Pr"]["group"]     = 3;
        elements["Pr"]["matType"]   = 9;
        elements["Pr"]["atomNum"]   = 59;
        elements["Pr"]["val"]       = 4;
        //elements["Po"]["vdw"]       = 1.85;
        //elements["Po"]["cova"]      = 1.60;
        elemProps["Pr"]["vdw"]      = 1.82;  //atom radius
        elemProps["Pr"]["cova"]     = 2.03;
        elemProps["Pr"]["ionM+"]    = 1.08;


        elements["Nd"]["row"]       = 6;
        elements["Nd"]["group"]     = 3;
        elements["Nd"]["matType"]   = 9;
        elements["Nd"]["atomNum"]   = 60;
        elements["Nd"]["val"]       = 3;
        //elements["Po"]["vdw"]       = 1.85;
        //elements["Po"]["cova"]      = 1.60;
        elemProps["Nd"]["vdw"]      = 1.81;  //atom radius
        elemProps["Nd"]["cova"]     = 2.01;
        elemProps["Nd"]["ionM+"]    = 1.08;


        elements["Pm"]["row"]       = 6;
        elements["Pm"]["group"]     = 3;
        elements["Pm"]["matType"]   = 9;
        elements["Pm"]["atomNum"]   = 61;
        elements["Pm"]["val"]       = 2;
        //elements["Po"]["vdw"]       = 1.85;
        //elements["Po"]["cova"]      = 1.60;
        elemProps["Pm"]["vdw"]      = 1.59;
        elemProps["Pm"]["cova"]     = 1.87;
        elemProps["Pm"]["ionM+"]    = 1.08;


        elements["Sm"]["row"]       = 6;
        elements["Sm"]["group"]     = 3;
        elements["Sm"]["matType"]   = 9;
        elements["Sm"]["atomNum"]   = 62;
        elements["Sm"]["val"]       = 3;
        //elements["Po"]["vdw"]       = 1.85;
        //elements["Po"]["cova"]      = 1.60;
        elemProps["Sm"]["vdw"]      = 1.59;
        elemProps["Sm"]["cova"]     = 1.87;
        elemProps["Sm"]["ionM+"]    = 1.08;

        elements["Eu"]["row"]       = 6;
        elements["Eu"]["group"]     = 3;
        elements["Eu"]["matType"]   = 9;
        elements["Eu"]["atomNum"]   = 63;
        elements["Eu"]["val"]       = 6;
        //elements["Po"]["vdw"]       = 1.85;
        //elements["Po"]["cova"]      = 1.60;
        elemProps["Eu"]["vdw"]      = 1.59;
        elemProps["Eu"]["cova"]     = 1.87;
        elemProps["Eu"]["ionM+"]    = 1.08;

        elements["Gd"]["row"]       = 6;
        elements["Gd"]["group"]     = 3;
        elements["Gd"]["matType"]   = 9;
        elements["Gd"]["atomNum"]   = 64;
        elements["Gd"]["val"]       = 3;
        elemProps["Gd"]["vdw"]      = 1.80;
        elemProps["Gd"]["cova"]     = 1.96;
        elemProps["Gd"]["ionM+"]    = 1.08;

        elements["Tb"]["row"]       = 6;
        elements["Tb"]["group"]     = 3;
        elements["Tb"]["matType"]   = 9;
        elements["Tb"]["atomNum"]   = 65;
        elements["Tb"]["val"]       = 4;
        elemProps["Tb"]["vdw"]      = 1.77;
        elemProps["Tb"]["cova"]     = 1.94;
        elemProps["Tb"]["ionM+"]    = 1.08;

        elements["Dy"]["row"]       = 6;
        elements["Dy"]["group"]     = 3;
        elements["Dy"]["matType"]   = 9;
        elements["Dy"]["atomNum"]   = 66;
        elements["Dy"]["val"]       = 3;
        elemProps["Dy"]["vdw"]      = 1.78;
        elemProps["Dy"]["cova"]     = 1.94;
        elemProps["Dy"]["ionM+"]    = 1.08;

        elements["Ho"]["row"]       = 6;
        elements["Ho"]["group"]     = 3;
        elements["Ho"]["matType"]   = 9;
        elements["Ho"]["atomNum"]   = 67;
        elements["Ho"]["val"]       = 3;
        elemProps["Ho"]["vdw"]      = 1.76;
        elemProps["Ho"]["cova"]     = 1.92;
        elemProps["Ho"]["ionM+"]    = 1.08;

        elements["Er"]["row"]       = 6;
        elements["Er"]["group"]     = 3;
        elements["Er"]["matType"]   = 9;
        elements["Er"]["atomNum"]   = 68;
        elements["Er"]["val"]       = 3;
        elemProps["Er"]["vdw"]      = 1.76;
        elemProps["Er"]["cova"]     = 1.89;
        elemProps["Er"]["ionM+"]    = 1.08;

        elements["Tm"]["row"]       = 6;
        elements["Tm"]["group"]     = 3;
        elements["Tm"]["matType"]   = 9;
        elements["Tm"]["atomNum"]   = 69;
        elements["Tm"]["val"]       = 3;
        elemProps["Tm"]["vdw"]      = 1.76;
        elemProps["Tm"]["cova"]     = 1.90;
        elemProps["Tm"]["ionM+"]    = 1.08;

        elements["Yb"]["row"]       = 6;
        elements["Yb"]["group"]     = 3;
        elements["Yb"]["matType"]   = 9;
        elements["Yb"]["atomNum"]   = 70;
        elements["Yb"]["val"]       = 1;
        elemProps["Yb"]["vdw"]      = 1.76;
        elemProps["Yb"]["cova"]     = 1.87;
        elemProps["Yb"]["ionM+"]    = 1.08;

        elements["Lu"]["row"]       = 6;
        elements["Lu"]["group"]     = 3;
        elements["Lu"]["matType"]   = 9;
        elements["Lu"]["atomNum"]   = 71;
        elements["Lu"]["val"]       = 1;
        elemProps["Lu"]["vdw"]      = 1.74;
        elemProps["Lu"]["cova"]     = 1.87;
        elemProps["Lu"]["ionM+"]    = 1.08;

        elements["Ac"]["row"]       = 7;
        elements["Ac"]["group"]     = 3;
        elements["Ac"]["matType"]   = 9;
        elements["Ac"]["atomNum"]   = 89;
        elements["Ac"]["val"]       = 3;
        elemProps["Ac"]["vdw"]      = 1.76;
        elemProps["Ac"]["cova"]     = 2.15;
        elemProps["Ac"]["ionM+"]    = 1.08;

        elements["Th"]["row"]       = 7;
        elements["Th"]["group"]     = 3;
        elements["Th"]["matType"]   = 9;
        elements["Th"]["atomNum"]   = 90;
        elements["Th"]["val"]       = 4;
        elemProps["Th"]["vdw"]      = 1.80;
        elemProps["Th"]["cova"]     = 2.06;
        elemProps["Th"]["ionM+"]    = 1.08;

        elements["Pa"]["row"]       = 7;
        elements["Pa"]["group"]     = 3;
        elements["Pa"]["matType"]   = 9;
        elements["Pa"]["atomNum"]   = 91;
        elements["Pa"]["val"]       = 4;
        elemProps["Pa"]["vdw"]      = 1.63;
        elemProps["Pa"]["cova"]     = 2.00;
        elemProps["Pa"]["ionM+"]    = 1.08;

        elements["U"]["row"]       = 7;
        elements["U"]["group"]     = 3;
        elements["U"]["matType"]   = 9;
        elements["U"]["atomNum"]   = 92;
        elements["U"]["val"]       = 6;
        elemProps["U"]["vdw"]      = 1.86;
        elemProps["U"]["cova"]     = 1.96;
        elemProps["U"]["ionM+"]    = 1.08;

        elements["Np"]["row"]       = 7;
        elements["Np"]["group"]     = 3;
        elements["Np"]["matType"]   = 9;
        elements["Np"]["atomNum"]   = 93;
        elements["Np"]["val"]       = 4;
        elemProps["Np"]["vdw"]      = 1.55;
        elemProps["Np"]["cova"]     = 1.90;
        elemProps["Np"]["ionM+"]    = 1.08;

        elements["Pu"]["row"]       = 7;
        elements["Pu"]["group"]     = 3;
        elements["Pu"]["matType"]   = 9;
        elements["Pu"]["atomNum"]   = 94;
        elements["Pu"]["val"]       = 6;
        elemProps["Pu"]["vdw"]      = 1.59;
        elemProps["Pu"]["cova"]     = 1.87;
        elemProps["Pu"]["ionM+"]    = 1.08;

        elements["Am"]["row"]       = 7;
        elements["Am"]["group"]     = 3;
        elements["Am"]["matType"]   = 9;
        elements["Am"]["atomNum"]   = 95;
        elements["Am"]["val"]       = 6;
        elemProps["Am"]["vdw"]      = 1.73;
        elemProps["Am"]["cova"]     = 1.80;
        elemProps["Am"]["ionM+"]    = 1.08;

        elements["Cm"]["row"]       = 7;
        elements["Cm"]["group"]     = 3;
        elements["Cm"]["matType"]   = 9;
        elements["Cm"]["atomNum"]   = 96;
        elements["Cm"]["val"]       = 4;
        elemProps["Cm"]["vdw"]      = 1.74;
        elemProps["Cm"]["cova"]     = 1.69;
        elemProps["Cm"]["ionM+"]    = 1.08;


        elements["Bk"]["row"]       = 7;
        elements["Bk"]["group"]     = 3;
        elements["Bk"]["matType"]   = 9;
        elements["Bk"]["atomNum"]   = 97;
        elements["Bk"]["val"]       = 4;
        elemProps["Bk"]["vdw"]      = 1.70;
        elemProps["Bk"]["cova"]     = 1.80;
        elemProps["Bk"]["ionM+"]    = 1.08;

        // The following are  rare-earth synthetic elements
        // where atom, vdw and covalent radii are not available
        // Put those props here just as position holders

        elements["Cf"]["row"]       = 7;
        elements["Cf"]["group"]     = 3;
        elements["Cf"]["matType"]   = 9;
        elements["Cf"]["atomNum"]   = 98;
        elements["Cf"]["val"]       = 4;
        elemProps["Cf"]["vdw"]      = 1.73;       // not sure
        elemProps["Cf"]["cova"]     = 1.80;       // not sure
        elemProps["Cf"]["ionM+"]    = 1.08;


        elements["Es"]["row"]       = 7;
        elements["Es"]["group"]     = 3;
        elements["Es"]["matType"]   = 9;
        elements["Es"]["atomNum"]   = 99;
        elements["Es"]["val"]       = 3;
        elemProps["Es"]["vdw"]      = 1.73;       // not sure
        elemProps["Es"]["cova"]     = 1.80;       // not sure
        elemProps["Es"]["ionM+"]    = 1.08;


        elements["Fm"]["row"]       = 7;
        elements["Fm"]["group"]     = 3;
        elements["Fm"]["matType"]   = 9;
        elements["Fm"]["atomNum"]   = 100;
        elements["Fm"]["val"]       = 3;
        elemProps["Fm"]["vdw"]      = 1.73;       // not sure
        elemProps["Fm"]["cova"]     = 1.80;       // not sure
        elemProps["Fm"]["ionM+"]    = 1.08;

        elements["Md"]["row"]       = 7;
        elements["Md"]["group"]     = 3;
        elements["Md"]["matType"]   = 9;
        elements["Md"]["atomNum"]   = 101;
        elements["Md"]["val"]       = 3;
        elemProps["Md"]["vdw"]      = 1.73;       // not sure
        elemProps["Md"]["cova"]     = 1.80;       // not sure
        elemProps["Md"]["ionM+"]    = 1.08;

        elements["No"]["row"]       = 7;
        elements["No"]["group"]     = 3;
        elements["No"]["matType"]   = 9;
        elements["No"]["atomNum"]   = 102;
        elements["No"]["val"]       = 3;
        elemProps["No"]["vdw"]      = 1.73;       // not sure
        elemProps["No"]["cova"]     = 1.80;       // not sure
        elemProps["No"]["ionM+"]    = 1.08;

        elements["Lr"]["row"]       = 7;
        elements["Lr"]["group"]     = 3;
        elements["Lr"]["matType"]   = 9;
        elements["Lr"]["atomNum"]   = 103;
        elements["Lr"]["val"]       = 3;
        elemProps["Lr"]["vdw"]      = 1.73;       // not sure
        elemProps["Lr"]["cova"]     = 1.80;       // not sure
        elemProps["Lr"]["ionM+"]    = 1.08;

        //some elements have multiple valences

        extraValences["P"].push_back(3);

        extraValences["N"].push_back(5);

        extraValences["S"].push_back(2);
        extraValences["S"].push_back(4);

        extraValences["Se"].push_back(2);
        extraValences["Se"].push_back(4);

        extraValences["Hg"].push_back(1);
        extraValences["Cu"].push_back(1);

        extraValences["Br"].push_back(3);
        extraValences["Br"].push_back(5);
        extraValences["Br"].push_back(7);

        for (std::map<ID, std::map<ID, int> >::iterator iEl=elements.begin();
                iEl!=elements.end(); iEl++)
        {
            if (iEl->second["group"]==1)
            {
                iEl->second["oxiNum"]= 1;
            }
            else if (iEl->second["group"]==2)
            {
                iEl->second["oxiNum"]=2;
            }
            else if (iEl->first == "O")
            {
                iEl->second["oxiNum"] = 2;
            }
            else if (iEl->first =="F")
            {
                iEl->second["oxiNum"]=1;
            }
            else
            {
                // Oxidation number of the rest of elements
                // are decided on-fly
                iEl->second["oxiNum"]=0;
            }

            // lone-pair electrons

            if (iEl->first == "N")
            {
                iEl->second["pairedVE"] = 2;
            }
            else if (iEl->first=="S")
            {
                iEl->second["pairedVE"] = 4;

            }
            else if (iEl->first=="P")
            {
                iEl->second["pairedVE"] = 2;
            }
            else if (iEl->first=="Cl" || iEl->first=="Br"
                    || iEl->first=="I" || iEl->first=="At")
            {
                iEl->second["pairedVE"] = 6;
            }
            else
            {
                // for metal and C atoms, set number of default lone-pair
                // electrons to zero

                iEl->second["pairedVE"] = 0;
            }

        }


        // Hydrogen
        electronConf["H"].push_back(std::make_pair ("1s",1));

        // Organic set or Non-metal
        electronConf["C"].push_back(std::make_pair ("2s",2));
        electronConf["C"].push_back(std::make_pair ("2p",2));

        electronConf["N"].push_back(std::make_pair ("2s",2));
        electronConf["N"].push_back(std::make_pair ("2p",3));

        electronConf["O"].push_back(std::make_pair ("2s",2));
        electronConf["O"].push_back(std::make_pair ("2p",4));


        electronConf["P"].push_back(std::make_pair ("3s",2));
        electronConf["P"].push_back(std::make_pair ("3p",3));

        electronConf["S"].push_back(std::make_pair ("3s",2));
        electronConf["S"].push_back(std::make_pair ("3p",4));


        electronConf["Se"].push_back(std::make_pair ("4s",2));
        electronConf["Se"].push_back(std::make_pair ("3d",10));
        electronConf["Se"].push_back(std::make_pair ("4p",4));

        // Halogens
        electronConf["F"].push_back(std::make_pair ("2s",2));
        electronConf["F"].push_back(std::make_pair ("2p",5));

        electronConf["Cl"].push_back(std::make_pair ("3s",2));
        electronConf["Cl"].push_back(std::make_pair ("3p",5));

        electronConf["Br"].push_back(std::make_pair ("4s",2));
        electronConf["Br"].push_back(std::make_pair ("4p",5));


        electronConf["I"].push_back(std::make_pair ("5s",2));
        electronConf["I"].push_back(std::make_pair ("5p",5));


        electronConf["At"].push_back(std::make_pair ("6s",2));
        electronConf["At"].push_back(std::make_pair ("6p",5));


        // Alkali-Metals
        electronConf["Li"].push_back(std::make_pair ("2s",1));

        electronConf["Na"].push_back(std::make_pair ("3s",1));

        electronConf["K"].push_back(std::make_pair ("4s",1));

        electronConf["Rb"].push_back(std::make_pair ("5s",1));

        electronConf["Cs"].push_back(std::make_pair ("6s",1));

        electronConf["Fr"].push_back(std::make_pair ("7s",1));

        // Alkaline-earth-Metals
        electronConf["Be"].push_back(std::make_pair ("2s",2));

        electronConf["Mg"].push_back(std::make_pair ("3s",2));

        electronConf["Ca"].push_back(std::make_pair ("4s",2));

        electronConf["Sr"].push_back(std::make_pair ("5s",2));

        electronConf["Ba"].push_back(std::make_pair ("6s",2));

        electronConf["Ra"].push_back(std::make_pair ("7s",2));

        // Transition-Metals

        electronConf["Sc"].push_back(std::make_pair ("3d",1));
        electronConf["Y"].push_back(std::make_pair ("4d",1));

        electronConf["Ti"].push_back(std::make_pair ("3d",2));
        electronConf["Zr"].push_back(std::make_pair ("4d",2));
        electronConf["Hf"].push_back(std::make_pair ("5d",2));
        electronConf["Rf"].push_back(std::make_pair ("6d",2));


        electronConf["V"].push_back(std::make_pair ("3d",3));


        electronConf["Nb"].push_back(std::make_pair ("5s",1));
        electronConf["Nb"].push_back(std::make_pair ("4d",4));
        electronConf["Ta"].push_back(std::make_pair ("5d",3));
        electronConf["Db"].push_back(std::make_pair ("6d",3));

        electronConf["Cr"].push_back(std::make_pair ("4s",1));
        electronConf["Cr"].push_back(std::make_pair ("3d",5));
        electronConf["Mo"].push_back(std::make_pair ("5s",1));
        electronConf["Mo"].push_back(std::make_pair ("4d",5));

        electronConf["W"].push_back(std::make_pair ("5d",4));
        electronConf["Sg"].push_back(std::make_pair ("6d",4));

        electronConf["Mn"].push_back(std::make_pair ("3d",5));
        electronConf["Tc"].push_back(std::make_pair ("4d",5));
        electronConf["Re"].push_back(std::make_pair ("5d",5));
        electronConf["Bh"].push_back(std::make_pair ("6d",5));

        electronConf["Fe"].push_back(std::make_pair ("3d",6));
        electronConf["Fe"].push_back(std::make_pair ("4s",2));
        electronConf["Ru"].push_back(std::make_pair ("5s",1));
        electronConf["Ru"].push_back(std::make_pair ("4d",7));
        electronConf["Os"].push_back(std::make_pair ("5d",6));
        electronConf["Hs"].push_back(std::make_pair ("6d",6));

        electronConf["Co"].push_back(std::make_pair ("3d",7));
        electronConf["Rh"].push_back(std::make_pair ("5s",1));
        electronConf["Rh"].push_back(std::make_pair ("4d",8));
        electronConf["Ir"].push_back(std::make_pair ("5d",7));
        electronConf["Mt"].push_back(std::make_pair ("6d",7));

        electronConf["Ni"].push_back(std::make_pair ("3d",8));
        electronConf["Pd"].push_back(std::make_pair ("4d",10));

        electronConf["Pt"].push_back(std::make_pair ("6s",1));
        electronConf["Pt"].push_back(std::make_pair ("5d",9));
        electronConf["Ds"].push_back(std::make_pair ("7s",1));
        electronConf["Ds"].push_back(std::make_pair ("6d",9));

        electronConf["Cu"].push_back(std::make_pair ("4s",1));
        electronConf["Cu"].push_back(std::make_pair ("3d",10));
        electronConf["Ag"].push_back(std::make_pair ("5s",1));
        electronConf["Ag"].push_back(std::make_pair ("4d",10));
        electronConf["Au"].push_back(std::make_pair ("6s",1));
        electronConf["Au"].push_back(std::make_pair ("5d",10));
        electronConf["Rg"].push_back(std::make_pair ("7s",1));
        electronConf["Rg"].push_back(std::make_pair ("6d",10));

        electronConf["Zn"].push_back(std::make_pair ("4s",2));
        electronConf["Zn"].push_back(std::make_pair ("3d",10));
        electronConf["Cd"].push_back(std::make_pair ("5s",2));
        electronConf["Cd"].push_back(std::make_pair ("4d",10));
        electronConf["Hg"].push_back(std::make_pair ("6s",2));
        electronConf["Hg"].push_back(std::make_pair ("5d",10));
        electronConf["Cn"].push_back(std::make_pair ("7s",2));
        electronConf["Cn"].push_back(std::make_pair ("6d",10));

        // Other-Metal
        electronConf["Al"].push_back(std::make_pair ("3s",2));
        electronConf["Al"].push_back(std::make_pair ("3p",1));
        electronConf["Ga"].push_back(std::make_pair ("4s",2));
        electronConf["Ga"].push_back(std::make_pair ("4p",1));
        electronConf["In"].push_back(std::make_pair ("5s",2));
        electronConf["In"].push_back(std::make_pair ("4d",10));
        electronConf["In"].push_back(std::make_pair ("5p",1));
        electronConf["Ti"].push_back(std::make_pair ("4f",14));
        electronConf["Ti"].push_back(std::make_pair ("5d",10));
        electronConf["Ti"].push_back(std::make_pair ("6p",1));
        electronConf["Sn"].push_back(std::make_pair ("5s",2));
        electronConf["Sn"].push_back(std::make_pair ("4d",10));
        electronConf["Sn"].push_back(std::make_pair ("5p",2));

        electronConf["Pb"].push_back(std::make_pair ("6s",2));
        electronConf["Pb"].push_back(std::make_pair ("4f",14));
        electronConf["Pb"].push_back(std::make_pair ("5d",10));
        electronConf["Pb"].push_back(std::make_pair ("6p",2));

        electronConf["Bi"].push_back(std::make_pair ("6s",2));
        electronConf["Bi"].push_back(std::make_pair ("4f",14));
        electronConf["Bi"].push_back(std::make_pair ("5d",10));
        electronConf["Bi"].push_back(std::make_pair ("6p",3));

        // Semimetallics

        electronConf["B"].push_back(std::make_pair ("2s",2));
        electronConf["B"].push_back(std::make_pair ("2p",1));
        electronConf["Si"].push_back(std::make_pair ("3s",2));
        electronConf["Si"].push_back(std::make_pair ("3p",2));
        electronConf["Ge"].push_back(std::make_pair ("4s",2));
        electronConf["Ge"].push_back(std::make_pair ("3d",10));
        electronConf["Ge"].push_back(std::make_pair ("4p",2));
        electronConf["As"].push_back(std::make_pair ("4s",2));
        electronConf["As"].push_back(std::make_pair ("3d",10));
        electronConf["As"].push_back(std::make_pair ("4p",3));
        electronConf["Sb"].push_back(std::make_pair ("5s",2));
        electronConf["Sb"].push_back(std::make_pair ("4d",10));
        electronConf["Sb"].push_back(std::make_pair ("5p",3));
        electronConf["Te"].push_back(std::make_pair ("5s",2));
        electronConf["Te"].push_back(std::make_pair ("4d",10));
        electronConf["Te"].push_back(std::make_pair ("5p",4));
        electronConf["Po"].push_back(std::make_pair ("6s",2));
        electronConf["Po"].push_back(std::make_pair ("4f",14));
        electronConf["Po"].push_back(std::make_pair ("5d",10));
        electronConf["Po"].push_back(std::make_pair ("6p",4));

        // Electronegativity of atoms, Allen scale
        // Allen, L.C., J. Am. Chem. Soc., 111, 9003, 1989
        // Unknown = -1.0

        // Hydrogen
        electronegativities["H"]  = 2.300;


        // Organic set or Non-metal
        electronegativities["C"]  = 2.544;
        electronegativities["N"]  = 3.066;
        electronegativities["O"]  = 3.610;
        electronegativities["P"]  = 2.253;
        electronegativities["S"]  = 2.589;
        electronegativities["Se"] = 2.434;

        // Halogens
        electronegativities["F"]  = 4.193;
        electronegativities["Cl"] = 2.869;
        electronegativities["Br"] = 2.685;
        electronegativities["I"]  = 2.359;
        electronegativities["At"] = 2.39;

        // Alkali-Metals
        electronegativities["Li"] = 0.912;
        electronegativities["Na"] = 0.869;
        electronegativities["K"]  = 0.734;
        electronegativities["Rb"] = 0.706;
        electronegativities["Cs"] = 0.659;
        electronegativities["Fr"] = 0.67;

        // Alkaline-earth-Metals
        electronegativities["Be"] = 1.576;
        electronegativities["Mg"] = 1.293;
        electronegativities["Ca"] = 1.034;
        electronegativities["Sr"] = 0.963;
        electronegativities["Ba"] = 0.881;
        electronegativities["Ra"] = 0.89;

        // Transition-Metals

        electronegativities["Sc"] = 1.19;
        electronegativities["Y"]  = 1.12;
        electronegativities["Ti"] = 1.38;
        electronegativities["Zr"] = 1.32;
        electronegativities["Hf"] = 1.16;
        electronegativities["Rf"] = -1.0;
        electronegativities["V"]  = 1.53;
        electronegativities["Nb"] = 1.41;
        electronegativities["Ta"] = 1.34;
        electronegativities["Db"] = -1.0;
        electronegativities["Cr"] = 1.65;
        electronegativities["Mo"] = 1.47;
        electronegativities["W"]  = 1.47;
        electronegativities["Sg"] = -1.0;
        electronegativities["Mn"] = 1.75;
        electronegativities["Tc"] = 1.51;
        electronegativities["Re"] = 1.60;
        electronegativities["Bh"] = -1.0;
        electronegativities["Fe"] = 1.80;
        electronegativities["Ru"] = 1.54;
        electronegativities["Os"] = 1.65;
        electronegativities["Hs"] = -1.0;
        electronegativities["Co"] = 1.84;
        electronegativities["Rh"] = 1.56;
        electronegativities["Ir"] = 1.68;
        electronegativities["Mt"] = -1.0;
        electronegativities["Ni"] = 1.88;
        electronegativities["Pd"] = 1.59;
        electronegativities["Pt"] = 1.72;
        electronegativities["Ds"] = -1.0;
        electronegativities["Cu"] = 1.85;
        electronegativities["Ag"] = 1.87;
        electronegativities["Au"] = 1.92;
        electronegativities["Rg"] = -1.0;
        electronegativities["Zn"] = 1.59;
        electronegativities["Cd"] = 1.52;
        electronegativities["Hg"] = 1.76;
        electronegativities["Cn"] = -1.0;

        // Other-Metal
        electronegativities["Al"] = 1.613;
        electronegativities["Ga"] = 1.756;
        electronegativities["In"] = 1.656;
        electronegativities["Ti"] = 1.789;
        electronegativities["Sn"] = 1.824;
        electronegativities["Pb"] = 1.854;
        electronegativities["Bi"] = 2.01;

        // Semimetallics

        electronegativities["B"]  = 2.051;
        electronegativities["Si"] = 1.916;
        electronegativities["Ge"] = 1.994;
        electronegativities["As"] = 2.211;
        electronegativities["Sb"] = 1.984;
        electronegativities["Te"] = 2.158;
        electronegativities["Po"] = 2.190;

    }

    PeriodicTable::~PeriodicTable()
    {
    }

    void PeriodicTable::compareIdealDists(std::map<std::string, std::map<
                                          std::string, std::map<std::string,
                                          REAL> > >  & tAllowedRangeDists,
                                          double tDelta)
    {
        // This function is for metal related bonds

        //std::ofstream compTab("idealComp.table");
        tDelta = tDelta/10.0;
        std::vector<ID>  metTab, metLoidTab;

        initMetalTab(metTab);
        initMetalloidTab(metLoidTab);

        for (std::map<ID,  std::map<ID, REAL> >::iterator
             iProp1=elemProps.begin(); iProp1 != elemProps.end(); iProp1++)
        {
            ID elem1 = iProp1->first;
            if (isMetal(metTab, elem1) || isMetalloid(metLoidTab, elem1))
            {
                for (std::map<ID,  std::map<ID, REAL> >::iterator
                     iProp2=elemProps.begin(); iProp2 != elemProps.end();
                     iProp2++)
                {

                    double d1=0.0, d2=0.0, dI=0;
                    ID elem2 = iProp2->first;
                    if (isMetal(metTab, elem2))
                    {
                        if (elemProps[elem2].find("ionM-")
                            !=elemProps[elem2].end())
                        {
                            dI = elemProps[elem1]["ionM+"]
                               + elemProps[elem2]["ionM-"];
                        }
                        else
                        {
                            dI = elemProps[elem1]["ionM+"]
                               + elemProps[elem2]["ionM+"];
                        }
                    }
                    else
                    {
                        if (elemProps[elem1].find("cova")
                            !=elemProps[elem1].end()
                            &&
                            elemProps[elem2].find("cova")
                            !=elemProps[elem2].end())
                        {
                            d1 = elemProps[elem1]["cova"]
                               + elemProps[elem2]["cova"];
                        }

                        if (elemProps[elem1].find("ionM+")
                             !=elemProps[elem1].end()
                             &&
                             elemProps[elem2].find("ionM-")
                             !=elemProps[elem2].end())
                        {
                            d2 = elemProps[elem1]["ionM+"]
                                + elemProps[elem2]["ionM-"];
                        }
                        /*
                        if(elem1.find("Mo")!=std::string::npos
                           && elem2.find("N")!=std::string::npos)
                        {
                            std::cout << "elem1 " << elem1 << " " << elem2 << std::endl;
                            std::cout << " d1 " << d1 << std::endl;
                            std::cout << " d2 " << d2 << std::endl;

                        }
                        */
                        if (d1 >= d2)
                        {
                            dI = d1;
                        }
                        else
                        {
                            dI = d2;
                        }
                    }

                    if (dI !=0)
                    {
                        tAllowedRangeDists[elem1][elem2]["max"]= dI*(1+tDelta);
                        tAllowedRangeDists[elem1][elem2]["min"]= dI*(1-tDelta);
                    }
                    else
                    {
                        std::cout << "The distance range between elements "
                                  << elem1 << " and " << elem2
                                  << " is not determined " << std::endl;
                    }
                }
            }
        }

        /*
        for (std::map<std::string, std::map<std::string,
             std::map<std::string, REAL> > >::iterator
             iR1=tAllowedRangeDists.begin();
             iR1 != tAllowedRangeDists.end(); iR1++)
        {
            for (std::map<std::string, std::map<std::string, REAL> >::iterator
                 iR2=iR1->second.begin(); iR2 !=iR1->second.end(); iR2++)
            {
                std::cout << iR1->first << "\t"
                          << iR2->first << "\t"
                          << iR2->second["min"] << "\t"
                          << iR2->second["max"] << "\n";
            }
        }
        */



    }

    void PeriodicTable::construAtomIdealDists(std::map<std::string,
       std::map<std::string,std::map<std::string,REAL> > >& tAllowedRangeDists,
       double tDelta)
    {
        const double dS = 1.0, dL=2.0, dDiff=0.5;

        // Try use electro-negative to decide the ideal pair distances.
        for (std::map<ID,  std::map<ID, REAL> >::iterator
             iProp1=elemProps.begin(); iProp1 != elemProps.end(); iProp1++)
        {
            ID elem1 = iProp1->first;
            for (std::map<ID,  std::map<ID, REAL> >::iterator
                     iProp2=elemProps.begin(); iProp2 != elemProps.end();
                     iProp2++)
            {
                ID elem2 = iProp2->first;
                double dIon, dCova, dMax, dIdeal=0.0;

                if (electronegativities.find(elem1) !=electronegativities.end()
                    &&
                    electronegativities.find(elem2) !=electronegativities.end())
                {
                    // dIon
                    if(electronegativities[elem1] <=
                       electronegativities[elem2])
                    {
                        if (elemProps[elem2].find("ionM-")
                            !=elemProps[elem2].end()
                            && elemProps[elem1].find("ionM+")
                            !=elemProps[elem1].end())
                        {
                            dIon=elemProps[elem1]["ionM+"]
                                +elemProps[elem2]["ionM-"];
                        }
                        else
                        {
                            dIon = elemProps[elem1]["cova"]
                                   +elemProps[elem2]["cova"];
                        }
                    }
                    else
                    {
                        if (elemProps[elem1].find("ionM-")
                            !=elemProps[elem1].end()
                            && elemProps[elem2].find("ionM+")
                            !=elemProps[elem2].end())
                        {
                            dIon=elemProps[elem1]["ionM-"]
                                +elemProps[elem2]["ionM+"];
                        }
                        else
                        {
                            dIon = elemProps[elem1]["cova"]
                                   +elemProps[elem2]["cova"];
                        }
                    }

                    // dCova
                    dCova = elemProps[elem1]["cova"]
                           +elemProps[elem2]["cova"];

                    // dMax
                    if (dIon >= dCova)
                    {
                        dMax = dIon;
                    }
                    else
                    {
                        dMax = dCova;
                    }


                    if (electronegativities[elem1] <= dS &&
                        electronegativities[elem2] <= dS)
                    {
                        dIdeal = dMax;
                    }
                    else if (electronegativities[elem1] >= dL &&
                             electronegativities[elem2] >= dL)
                    {
                        dIdeal = dCova;
                    }
                    else
                    {
                        double diffNE = fabs(electronegativities[elem1]
                                            -electronegativities[elem2]);
                        if (diffNE<= dDiff)
                        {
                            dIdeal = dCova;
                        }
                        else
                        {
                            dIdeal = dIon;
                        }
                    }

                    if (dIdeal !=0.0)
                    {
                        tAllowedRangeDists[elem1][elem2]["max"]
                                                  =dIdeal*(1+tDelta);
                        tAllowedRangeDists[elem1][elem2]["min"]
                                                  = dIdeal*(1-tDelta);
                    }
                    else
                    {
                        std::cout << "The distance range between elements "
                                  << elem1 << " and " << elem2
                                  << " is not determined " << std::endl;
                    }

                }
                else
                {
                    // dCova
                    dCova = elemProps[elem1]["cova"]
                           +elemProps[elem2]["cova"];


                }
            }
        }

    }

}
