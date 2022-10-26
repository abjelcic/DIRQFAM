!======================================================================!

      subroutine read_input()

!======================================================================!
      use fam_input;
      implicit none;
      integer           :: unit;
      integer           :: ios;
      character(len=32) :: nucName;
      character(len=32) :: parName;
      integer           :: nucA;


          open( status='old' , action='read' , form='formatted' , newunit=unit , file='dirqfam.dat' , iostat=ios );
          if( ios /= 0 ) stop 'Opening file dirqfam.dat failed.';

          read(unit,'(10x,2i5)')    n0f, n0b;
          read(unit,'(10x,f12.3)')  beta0;
          read(unit,*);
          read(unit,*);
          read(unit,'(a2,i4)')      nucName , nucA;
          read(unit,*);
          read(unit,*);
          read(unit,*);
          read(unit,'(10x,a10)')    parName;
          read(unit,*);
          read(unit,*);
          read(unit,*);
          read(unit,*);
          read(unit,*);
          read(unit,*);
          read(unit,*);
          read(unit,*);
          read(unit,'(19x,i12)')    calculation_type;
          read(unit,*);
          read(unit,*);
          read(unit,*);
          read(unit,*);
          read(unit,'(19x,i12)')    include_coulomb;
          read(unit,'(19x,i12)')    include_pairing;
          read(unit,'(19x,i12)')    NGH;
          read(unit,'(19x,i12)')    NGL;
          read(unit,'(19x,1f12.3)') gamma_smear;
          read(unit,'(19x,1f12.3)') selfConsistencyTolerance;
          read(unit,'(19x,i12)')    NoArnoldiVectors;
          read(unit,*);
          read(unit,'(19x,i12)')    J_multipole;
          read(unit,'(19x,i12)')    K_multipole;
          read(unit,'(19x,i12)')    Isospin;
          read(unit,*);
          read(unit,'(19x,1f12.3)') omega_start;
          read(unit,'(19x,1f12.3)') omega_end;
          read(unit,'(19x,1f12.3)') delta_omega;
          read(unit,*);
          read(unit,'(19x,1f12.3)') omega_print;
          read(unit,*);
          read(unit,'(19x,1f12.3)') omega_center;
          read(unit,'(19x,1f12.3)') omega_radius;
          read(unit,'(19x,i12)')    NoContourPoints;

          close(unit);


          nucleusName         = trim(adjustl(nucName));
          nucleusZ            = getAtomicNumberFromElementName( nucleusName );
          nucleusN            = nucA - nucleusZ;
          LagrangianModelName = trim(adjustl(parName));


          call assert( nucleusZ>=1 .and. nucleusN>=1 , 'Number of protons/neutrons wrong.' );

          call assert( all([mod(nucleusZ,2),mod(nucleusN,2)]==[0,0]) , 'Only even-even nuclei allowed.'   );
          call assert( all([mod(n0f     ,2),mod(n0b     ,2)]==[0,0]) , 'n0f and n0b must be even numbers.');

          call assert( LagrangianModelName=='DD-PC1' .or. LagrangianModelName=='DD-ME2' , 'Only DD-PC1 and DD-ME2 supported.' );


      return;

      contains

          integer function getAtomicNumberFromElementName( elementName ) result(Z)
              implicit none;
              character(len=*) , intent(in) :: elementName;

              Z = -1;

              if( elementName=='_H' ) Z =   1;
              if( elementName=='He' ) Z =   2;
              if( elementName=='Li' ) Z =   3;
              if( elementName=='Be' ) Z =   4;
              if( elementName=='_B' ) Z =   5;
              if( elementName=='_C' ) Z =   6;
              if( elementName=='_N' ) Z =   7;
              if( elementName=='_O' ) Z =   8;
              if( elementName=='_F' ) Z =   9;
              if( elementName=='Ne' ) Z =  10;
              if( elementName=='Na' ) Z =  11;
              if( elementName=='Mg' ) Z =  12;
              if( elementName=='Al' ) Z =  13;
              if( elementName=='Si' ) Z =  14;
              if( elementName=='_P' ) Z =  15;
              if( elementName=='_S' ) Z =  16;
              if( elementName=='Cl' ) Z =  17;
              if( elementName=='Ar' ) Z =  18;
              if( elementName=='_K' ) Z =  19;
              if( elementName=='Ca' ) Z =  20;
              if( elementName=='Sc' ) Z =  21;
              if( elementName=='Ti' ) Z =  22;
              if( elementName=='_V' ) Z =  23;
              if( elementName=='Cr' ) Z =  24;
              if( elementName=='Mn' ) Z =  25;
              if( elementName=='Fe' ) Z =  26;
              if( elementName=='Co' ) Z =  27;
              if( elementName=='Ni' ) Z =  28;
              if( elementName=='Cu' ) Z =  29;
              if( elementName=='Zn' ) Z =  30;
              if( elementName=='Ga' ) Z =  31;
              if( elementName=='Ge' ) Z =  32;
              if( elementName=='As' ) Z =  33;
              if( elementName=='Se' ) Z =  34;
              if( elementName=='Br' ) Z =  35;
              if( elementName=='Kr' ) Z =  36;
              if( elementName=='Rb' ) Z =  37;
              if( elementName=='Sr' ) Z =  38;
              if( elementName=='_Y' ) Z =  39;
              if( elementName=='Zr' ) Z =  40;
              if( elementName=='Nb' ) Z =  41;
              if( elementName=='Mo' ) Z =  42;
              if( elementName=='Tc' ) Z =  43;
              if( elementName=='Ru' ) Z =  44;
              if( elementName=='Rh' ) Z =  45;
              if( elementName=='Pd' ) Z =  46;
              if( elementName=='Ag' ) Z =  47;
              if( elementName=='Cd' ) Z =  48;
              if( elementName=='In' ) Z =  49;
              if( elementName=='Sn' ) Z =  50;
              if( elementName=='Sb' ) Z =  51;
              if( elementName=='Te' ) Z =  52;
              if( elementName=='_I' ) Z =  53;
              if( elementName=='Xe' ) Z =  54;
              if( elementName=='Cs' ) Z =  55;
              if( elementName=='Ba' ) Z =  56;
              if( elementName=='La' ) Z =  57;
              if( elementName=='Ce' ) Z =  58;
              if( elementName=='Pr' ) Z =  59;
              if( elementName=='Nd' ) Z =  60;
              if( elementName=='Pm' ) Z =  61;
              if( elementName=='Sm' ) Z =  62;
              if( elementName=='Eu' ) Z =  63;
              if( elementName=='Gd' ) Z =  64;
              if( elementName=='Tb' ) Z =  65;
              if( elementName=='Dy' ) Z =  66;
              if( elementName=='Ho' ) Z =  67;
              if( elementName=='Er' ) Z =  68;
              if( elementName=='Tm' ) Z =  69;
              if( elementName=='Yb' ) Z =  70;
              if( elementName=='Lu' ) Z =  71;
              if( elementName=='Hf' ) Z =  72;
              if( elementName=='Ta' ) Z =  73;
              if( elementName=='_W' ) Z =  74;
              if( elementName=='Re' ) Z =  75;
              if( elementName=='Os' ) Z =  76;
              if( elementName=='Ir' ) Z =  77;
              if( elementName=='Pt' ) Z =  78;
              if( elementName=='Au' ) Z =  79;
              if( elementName=='Hg' ) Z =  80;
              if( elementName=='Tl' ) Z =  81;
              if( elementName=='Pb' ) Z =  82;
              if( elementName=='Bi' ) Z =  83;
              if( elementName=='Po' ) Z =  84;
              if( elementName=='At' ) Z =  85;
              if( elementName=='Rn' ) Z =  86;
              if( elementName=='Fr' ) Z =  87;
              if( elementName=='Ra' ) Z =  88;
              if( elementName=='Ac' ) Z =  89;
              if( elementName=='Th' ) Z =  90;
              if( elementName=='Pa' ) Z =  91;
              if( elementName=='_U' ) Z =  92;
              if( elementName=='Np' ) Z =  93;
              if( elementName=='Pu' ) Z =  94;
              if( elementName=='Am' ) Z =  95;
              if( elementName=='Cm' ) Z =  96;
              if( elementName=='Bk' ) Z =  97;
              if( elementName=='Cf' ) Z =  98;
              if( elementName=='Es' ) Z =  99;
              if( elementName=='Fm' ) Z = 100;
              if( elementName=='Md' ) Z = 101;
              if( elementName=='No' ) Z = 102;
              if( elementName=='Lr' ) Z = 103;
              if( elementName=='Rf' ) Z = 104;
              if( elementName=='Db' ) Z = 105;
              if( elementName=='Sg' ) Z = 106;
              if( elementName=='Bh' ) Z = 107;
              if( elementName=='Hs' ) Z = 108;
              if( elementName=='Mt' ) Z = 109;
              if( elementName=='Ds' ) Z = 110;
              if( elementName=='Rg' ) Z = 111;
              if( elementName=='Cn' ) Z = 112;
              if( elementName=='Nh' ) Z = 113;
              if( elementName=='Fl' ) Z = 114;
              if( elementName=='Mc' ) Z = 115;
              if( elementName=='Lv' ) Z = 116;
              if( elementName=='Ts' ) Z = 117;
              if( elementName=='Og' ) Z = 118;

              call assert( Z > 0 , 'Element name: ['//elementName//'] unknown.' );

              return;
          end function getAtomicNumberFromElementName

      end subroutine read_input
