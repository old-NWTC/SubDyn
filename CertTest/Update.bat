@ECHO OFF
@ECHO.


REM  Set up environment variables.  You will probably have to change these.


@SET ResultsDir=NREL_RESULTS

::=======================================================================================================


:CertTest


REM  SubDyn test sequence definition:




@CALL :CopySubDyn 01 out
@CALL :CopySubDyn 02 out
@CALL :CopySubDyn 03 out
@CALL :CopySubDyn 04 out
@CALL :CopySubDyn 05 out

goto END

rem ******************************************************
:CopySubDyn
:: Copy files to NREL_Results.

COPY Test%1\Test%1.SD.out Test%1\%ResultsDir%\Test%1.SD.out
COPY Test%1\Test%1.SD.sum Test%1\%ResultsDir%\Test%1.SD.sum


EXIT /B


:__ErrorExit
rem Creates a syntax error, stops immediately
()
EXIT /B


:END

EXIT /B