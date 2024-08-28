/*
 * File:   PDBFile.h
 * Author: flong
 *
 * Created on August 31, 2011, 11:44 AM
 */

#ifndef PDBFILE_H
#define	PDBFILE_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef PDBHEADERSEC_H
#include "PDBHeaderSec.h"
#endif

#ifndef FILE_H
#include "file.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef LIBG_ATOMASSEMBLY_H
#include "atomAssembly.h"
#endif

#ifndef RESIDUE_H
#include "residue.h"
#endif

#ifndef CHAIN_H
#include "chain.h"
#endif

#ifndef LIBG_MODEL_H
#include "libgmodel.h"
#endif

#ifndef BOND_H
#include "bond.h"
#endif


#ifndef ANGLE_H
#include "angle.h"
#endif

#ifndef SSBOND_H
#include "ssbond.h"
#endif

#ifndef LIBG_LINK_H
#include "libglink.h"
#endif

#ifndef SECONDARYSTRUCTURES_H
#include "secondaryStructures.h"
#endif

#ifndef MONOMERLIB_H
#include "MonomerLib.h"
#endif

#ifndef CRYSTINFO_H
#include "crystInfo.h"
#endif

#ifndef CHEMPROPSET_H
#include "chemPropSet.h"
#endif

#ifndef UTILITY_H
#include "utility.h"
#endif

namespace LIBMOL
{
    class Atom;
    class AtomDict;

    class Residue;
    class ModRes;
    class Chain;
    class Model;

    class Helix;
    class Sheet;
    class Strand;
    class Turn;

    class SSBond;
    class Link;

    class CrystInfo;

    enum PDBRecordType
    {
        UnknownType = 0,
        AnisouType,
        AtomType,
        AuthorType,
        CaveatType,
        CispepType,
        CompndType,
        ConectType,
        Cryst1Type,
        DbRefType,
        DbRef1Type,
        DbRef2Type,
        EndType,
        EndMdlType,
        ExpDtaType,
        ForMulType,
        FtNoteType,
        HeadType,
        HelixType,
        HetType,
        HetatmType,
        HetNamType,
        HetSynType,
        HydBndType,
        JrnlType,
        keywdsType,
        LinkType,
        MasterType,
        ModelType,
        ModResType,
        MTRIX1Type,
        MTRIX2Type,
        MTRIX3Type,
        ObslteType,
        ORIGX1Type,
        ORIGX2Type,
        ORIGX3Type,
        RemarkType,
        RevDatType,
        SCALE1Type,
        SCALE2Type,
        SCALE3Type,
        SeqAdvType,
        SeqResType,
        SheetType,
        SigAtmType,
        SigUijType,
        SiteType,
        SltBrgType,
        SourceType,
        SprsdeType,
        SsBondType,
        TerType,
        TitleType,
        TurnType,
        TVectType
    };


    class PDBFile : public File
    {
    public :

        // default constructor
        PDBFile();

        // copy constructor
        PDBFile(Name               tFname,
                std::ios_base::openmode tOpenMode);

        PDBFile(FileName            tFname,
                std::ios_base::openmode tOpenMode);

        // destructor
        virtual ~PDBFile();



        std::vector<std::string>  PDBRecordHeader;

        void setPDBRecordHeader();

        // structure that stores the records in the coordinate section.
        // Maybe using class Atom directly is simple enough
        struct PDBAtom
        {
           Name               recordName;
           SeriNumber         seriNum;
           std::string        name;
           std::string        altLoc;    // alternate location indicator
           std::string        resName;   // amino acid abbreviation for an atom
           ID                 atomID;        // atom role name. See PDB convention
           ID                 chainID;   // the chain ID for an atom
           SeriNumber         segNum;    // segment serial number
           std::string        insCode;   //  Code for insertion of residues
           std::string        elementType; // the element type of an atom

           std::string        itsSegID;       // segment identifier

        };


        // Extract information from a PDB file

        // Read line by line from a PDB file and
        // and get all information needed from it
        void setupSystem();

        // Get a record type from the first 6 letters in a PDB line
        int  getRecordType(ALine         tLine);

        /* The following member functions get different
           information from different type of record-lines
        */
        bool extractTitleInfo(ALine          tLine);
        // Primary Structure Section
        bool extractSEQRESInfo(ALine         tLine); // get chains
        bool extractMODRESInfo(ALine         tLine);

        // Secondary Structure Section
        bool extractHelixInfo(ALine          tLine);
        bool extractSheetInfo(ALine          tLine);

         // Connectivity Annotation Section
        bool extractSSBondInfo(ALine          tLine);

        bool extractLinkInfo(ALine          tLine);

        // Crystallographic and Coordinate Transformation Section
        bool extractCryst1Info(ALine         tLine);
        bool extractOrigXNInfo(ALine         tLine);
        bool extractScaleNInfo(ALine         tLine);
        bool extractMatrixNInfo(ALine         tLine);

        // Coordinate Section
        void initOneModel();
        void initOneModel(ALine               tLine);

        bool extractAtomInfo(ALine          tLine);
        bool extractAtomInfo(Atom      &    tAtom,
                             Residue   &    tResidue,
                             Chain     &    tChain,
                             ALine          tLine);

        bool extractAnisouInfo(ALine          tLine);
        bool extractTerInfo(ALine          tLine);
        bool extractHetatmInfo(ALine          tLine);
        bool extractEndMdlInfo(ALine          tLine);


        bool extractConnectInfo(ALine          tLine);

        bool extractBookeepInfo(ALine          tLine);

        // access some general information
        ID   getPDBID() const;

        Date getDepositedDate() const;

        std::string getDesciption();

        ObsoleteInfo getObsoleteInfo();

        std::vector<Model> getModelList();
        void               setModelList();

        void setResidueGroupType(FileName tF);

        // output a structure to a PDB file

        void writeTitleSection();
        void writePrimaryStructureSection();
        void writeHeterogenSection();
        void writeSecondaryStructureSection();
        void writeConnectivityAnnotationSection();
        void writeMiscellaneousFeaturesSection();
        void writeCrystallographicSection();
        void writeCoordinateSection();
        void writeConnectivitySection();
        void writeBookKeepingSection();

        // related bonds, bond angle, torsion angles, residues, chain. models
        std::vector<Model>      allModels;
        std::vector<Chain>      allChains;
        std::vector<ModRes>     allModRes;
        std::vector<Atom>       allAtomList;
        std::vector<Atom>       allHetAtmList;
        std::vector<SSBond>     allSSBonds;
        std::vector<Link>       allLinks;
        std::vector<Helix>      allHelices;
        std::vector<Sheet>      allSheets;

        HeadSec                 allHeaderInfo;

        std::ofstream           outFile;
        std::ifstream           inFile;

    private :


        PDBRecordType           itsTempRecordType;

        Model                *  itsCurrentModel;
        Atom                 *  itsCurrentAtom;
        bool                    itsCurrentAltLoc;
        Atom                 *  itsCurrentHetAtm;


        Residue              *  itsCurrentResidue;
        ModRes               *  itsCurrentModRes;
        Chain                *  itsCurrentChain;

        Helix                *  itsCurrentHelix;
        Sheet                *  itsCurrentSheet;


        SSBond               *  itsCurrentSSBond;

        CrystInfo            *   itsCrystInfo;

    };


    // A simple version of PDB processor, extract only atom information
    // in particular, atomic coordinates
    class DictPDBFile
    {
    public :

        // default constructor
        DictPDBFile();

        // copy constructor
        DictPDBFile(Name                    tFname,
                    std::ios_base::openmode tOpenMode);

        inline ~DictPDBFile()
        {
        }


        // Extract information from a PDB file

        // Read line by line from a PDB file and
        // and get all information needed from it
        void setupSystem();

        // atoms only
        std::map<ID, std::vector<REAL> >       allHetAtmList;

        std::ifstream                          inFile;


    };

    void extern outPDB(FileName tFName,
                       ID tMonoRootName,
                       std::vector<LIBMOL::AtomDict>& tAtoms,
                       int tMode);

}

#endif	/* PDBFILE_H */
