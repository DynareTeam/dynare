# Configuration file for building Dynare Windows Installer
# Uses "NullSoft Scriptable Installer System", aka NSIS (see http://nsis.sourceforge.net)
# NSIS can be run from both Windows and Linux (see "nsis" package in Debian)

# How to build the installer:
# - build: the preprocessor, the MEX binaries (for MATLAB and for Octave), and the documentation (PDF files + HTML manual)
# - run "makensis dynare.nsi" to create the installer
# - if there is no failure, this will create a file "dynare-VERSION-win.exe" in the current directory

!include dynare-version.nsi

SetCompressor /SOLID lzma

Name "Dynare ${VERSION}"

OutFile "dynare-${VERSION}-win.exe"

InstallDir "c:\dynare\${VERSION}"

# Use the Modern User Interface (version 2)
!include MUI2.nsh

!define MUI_WELCOMEPAGE_TEXT "This wizard will guide you through the installation of Dynare ${VERSION}.$\n$\nDynare is distributed under the GNU General Public License (GPL) version 3.$\n$\nIf you accept the license, click Next button to continue the installation."
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_COMPONENTS
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!define MUI_FINISHPAGE_LINK_LOCATION http://www.dynare.org
!define MUI_FINISHPAGE_LINK "Go to Dynare homepage"
!define MUI_FINISHPAGE_SHOWREADME $INSTDIR\README.txt
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_UNPAGE_WELCOME
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_UNPAGE_FINISH

!insertmacro MUI_LANGUAGE "English"

!define REGLOC "Software\Microsoft\Windows\CurrentVersion\Uninstall\Dynare ${VERSION}"
!define SMLOC "$SMPROGRAMS\Dynare ${VERSION}"

!macro DETERMINE_CONTEXT
 # Determine if we are admin or not
 # This will change the start menu directory and the registry root key (HKLM or HKLU)
 UserInfo::getAccountType
 Pop $0
 StrCmp $0 "Admin" +3
 SetShellVarContext current
 Goto +2
 SetShellVarContext all
!macroend

Section "Dynare core (preprocessor and M-files)"
 SectionIn RO
!insertmacro DETERMINE_CONTEXT
 SetOutPath $INSTDIR
 File README.txt ..\NEWS mexopts-win32.bat mexopts-win64.bat ..\license.txt

 SetOutPath $INSTDIR\matlab
 File /r ..\matlab\*.m

 SetOutPath $INSTDIR\matlab\preprocessor32
 File ..\matlab\preprocessor32\dynare_m.exe

 SetOutPath $INSTDIR\matlab\preprocessor64
 File ..\matlab\preprocessor64\dynare_m.exe

 SetOutPath $INSTDIR\contrib
 File /r ..\contrib\*.m

 SetOutPath $INSTDIR\scripts
 File /r ..\scripts\*

 WriteUninstaller $INSTDIR\uninstall.exe

 # Create start menu entries
 CreateDirectory "${SMLOC}"
 CreateShortcut "${SMLOC}\Uninstall.lnk" "$INSTDIR\uninstall.exe"

 # Create entry in "Add/Remove programs"
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "DisplayName" "Dynare ${VERSION}"
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "DisplayVersion" "${VERSION}"
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "InstallLocation" $INSTDIR
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "UninstallString" "$INSTDIR\uninstall.exe"
 WriteRegDWORD SHELL_CONTEXT "${REGLOC}" "NoModify" 1
 WriteRegDWORD SHELL_CONTEXT "${REGLOC}" "NoRepair" 1
SectionEnd

SectionGroup "MEX files for MATLAB"

Section "MEX files for MATLAB 32-bit, version 7.5 to 8.6 (R2007b to R2015b)"
 SetOutPath $INSTDIR\mex\matlab\win32-7.5-8.6
 File ..\mex\matlab\win32-7.5-8.6\*.mexw32
SectionEnd

Section "MEX files for MATLAB 64-bit, version 7.5 to 7.7 (R2007b to R2008b)"
 SetOutPath $INSTDIR\mex\matlab\win64-7.5-7.7
 File ..\mex\matlab\win64-7.5-7.7\*.mexw64
SectionEnd

Section "MEX files for MATLAB 64-bit, version 7.8 to 9.0 (R2009a to R2016a)"
 SetOutPath $INSTDIR\mex\matlab\win64-7.8-9.0
 File ..\mex\matlab\win64-7.8-9.0\*.mexw64
SectionEnd

SectionGroupEnd

Section "MEX files for Octave 3.6 (MinGW)"
 SetOutPath $INSTDIR\mex\octave
 File ..\mex\octave\*.mex ..\mex\octave\*.oct
SectionEnd

Section "Dynare++ (standalone executable)"
 SetOutPath $INSTDIR\dynare++
 File ..\dynare++\src\dynare++.exe ..\dynare++\extern\matlab\dynare_simul.m
SectionEnd

Section "Documentation and examples (Dynare and Dynare++)"
 SetOutPath $INSTDIR\doc
 File ..\doc\dynare.pdf ..\doc\guide.pdf ..\doc\userguide\UserGuide.pdf ..\doc\bvar-a-la-sims.pdf ..\doc\dr.pdf ..\doc\macroprocessor\macroprocessor.pdf ..\doc\preprocessor\preprocessor.pdf ..\doc\parallel\parallel.pdf ..\doc\gsa\gsa.pdf ..\doc\dseries-and-reporting\dseriesReporting.pdf

 SetOutPath $INSTDIR\doc\dynare.html
 File ..\doc\dynare.html\*.html ..\doc\dynare.html\*.png

 SetOutPath $INSTDIR\doc\dynare++
 File ..\dynare++\doc\dynare++-tutorial.pdf ..\dynare++\doc\dynare++-ramsey.pdf ..\dynare++\sylv\sylvester.pdf ..\dynare++\tl\cc\tl.pdf ..\dynare++\integ\cc\integ.pdf ..\dynare++\kord\kord.pdf

 CreateShortcut "${SMLOC}\Documentation.lnk" "$INSTDIR\doc"

 SetOutPath $INSTDIR\examples
 File ..\examples\*

 CreateShortcut "${SMLOC}\Examples.lnk" "$INSTDIR\examples"

SectionEnd

Section "Uninstall"
!insertmacro DETERMINE_CONTEXT

 # First delete the uninstaller
 Delete $INSTDIR\uninstall.exe
 Delete $INSTDIR\README.txt
 Delete $INSTDIR\NEWS
 Delete $INSTDIR\license.txt
 Delete $INSTDIR\mexopts-win32.bat
 Delete $INSTDIR\mexopts-win64.bat
 Rmdir /r $INSTDIR\matlab
 Rmdir /r $INSTDIR\contrib
 Rmdir /r $INSTDIR\mex
 Rmdir /r $INSTDIR\dynare++
 Rmdir /r $INSTDIR\doc
 Rmdir /r $INSTDIR\examples
 Rmdir /r $INSTDIR\scripts
 # We don't force deletion of installation directory (with /r), to avoid deleting important files
 Rmdir $INSTDIR

 # Delete start menu entries
 Rmdir /r "${SMLOC}"

 # Delete entry in "Add/Remove programs"
 DeleteRegKey SHELL_CONTEXT "${REGLOC}"
SectionEnd
