/* Author:  Bruce Swihart
   Date:  2013/02
   Warranty:  None

This .sas file accompanies the paper "Modeling sleep fragmentation in populations of sleep hypnograms".  */


* 5 STATE / 20 Transition Types:  Poisson format / loglinear modeling;



    
PROC IMPORT OUT= WORK.LLfull0
    DATAFILE="LLfull_allgroups_agesexracesmoking.csv"
    DBMS=CSV REPLACE;
    GETNAMES=YES;
    DATAROW=2;
RUN;



data LLFULL0;
    set LLFULL0;
    LOGTAR = -10;
    IF TAR GE 0.5 THEN LOGTAR = log(TAR);            /*Time at Risk (TAR) can only be in increments of 0.5 */
    IF TAR LT 0.5 THEN LOGTAR = log(0+.01);          /*So this adds a little to true 0s before log()*/  
    stagings235sum = staging2 + staging3 + staging5; /* sleep quality variable*/ 
run;

title "IND -- age sex race smoking";
PROC GENMOD data = LLfull0 ;

    class race smokstatus type group pptid / ref = FIRST param = REFERENCE ;
    model counts = type group type*group  age race sex smokstatus / D=POISSON LINK=LOG offset=LOGTAR ;

    estimate '0515_1R' group 1 0 0                                                                                                                               / exp ;
    estimate '0515_2R' group 1 0 0 type*group 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_SR' group 1 0 0 type*group 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_1W' group 1 0 0 type*group 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_2W' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_SW' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_R1' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_R2' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_RS' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_RW' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_W1' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_W2' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_WS' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_WR' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_S1' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_S2' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_12' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_1S' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0  / exp ;
    estimate '0515_21' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0  / exp ;
    estimate '0515_2S' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0  / exp ;
    
    estimate '1530_1R' group 0 1 0                                                                                                                               / exp ;
    estimate '1530_2R' group 0 1 0 type*group 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_SR' group 0 1 0 type*group 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_1W' group 0 1 0 type*group 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_2W' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_SW' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_R1' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_R2' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_RS' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_RW' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_W1' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_W2' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_WS' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_WR' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_S1' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_S2' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_12' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '1530_1S' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0  / exp ;
    estimate '1530_21' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0  / exp ;
    estimate '1530_2S' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0  / exp ;
    

    estimate '30++_1R' group 0 0 1                                                                                                                               / exp ;
    estimate '30++_2R' group 0 0 1 type*group 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_SR' group 0 0 1 type*group 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_1W' group 0 0 1 type*group 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_2W' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_SW' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_R1' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_R2' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_RS' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_RW' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_W1' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_W2' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_WS' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_WR' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_S1' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_S2' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_12' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0  / exp ;
    estimate '30++_1S' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0  / exp ;
    estimate '30++_21' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0  / exp ;
    estimate '30++_2S' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1  / exp ;
    

    REPEATED SUBJECT = PPTID / WITHINSUBJECT=TYPE CORR=IND SORTED MAXIT=10000;

RUN;




* 3 STATE /  6 Transition Types:  Poisson format / loglinear modeling;
PROC IMPORT OUT= WORK.LLfull
    DATAFILE="LLfull_allgroups_agesexracesmoking_0306.csv"
    DBMS=CSV REPLACE;
    GETNAMES=YES;
    DATAROW=2;
RUN;

title "IND -- 3/6  age sex race smoking";
PROC GENMOD data = LLfull ;

    class race smokstatus type group pptid / ref = FIRST param = REFERENCE ;
    model counts = type group type*group  age race sex smokstatus / D=POISSON LINK=LOG ;
    
    estimate '0515_NR' group 1 0 0                                          / exp ;
    estimate '0515_NW' group 1 0 0 type*group 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 / exp ;
    estimate '0515_RN' group 1 0 0 type*group 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 / exp ;
    estimate '0515_RW' group 1 0 0 type*group 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 / exp ;
    estimate '0515_WN' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 / exp ;
    estimate '0515_WR' group 1 0 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 / exp ;

    estimate '1530_NR' group 0 1 0                                          / exp ;
    estimate '1530_NW' group 0 1 0 type*group 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 / exp ;
    estimate '1530_RN' group 0 1 0 type*group 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 / exp ;
    estimate '1530_RW' group 0 1 0 type*group 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 / exp ;
    estimate '1530_WN' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 / exp ;
    estimate '1530_WR' group 0 1 0 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 / exp ;

    estimate '30++_NR' group 0 0 1                                          / exp ;
    estimate '30++_NW' group 0 0 1 type*group 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 / exp ;
    estimate '30++_RN' group 0 0 1 type*group 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 / exp ;
    estimate '30++_RW' group 0 0 1 type*group 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 / exp ;
    estimate '30++_WN' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 / exp ;
    estimate '30++_WR' group 0 0 1 type*group 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 / exp ;


    REPEATED SUBJECT = PPTID / WITHINSUBJECT=TYPE CORR=IND SORTED MAXIT=1000 CORRW;

RUN;


* 5 STATE / 20 Transition Types:  multistate survival;
PROC IMPORT OUT= WORK.SAfull
    DATAFILE= "SAfull_allgroups_agesexracesmoking.csv"
    DBMS=CSV REPLACE;
    GETNAMES=YES;
    DATAROW=2;
RUN;



title "SA 5/20";
PROC phreg data = SAfull covs(AGGREGATE) ;
    class race(REF=FIRST) smokstatus;
    model time*status(0) = 
        g1t01_1R g1t02_2R g1t03_SR
        g1t04_1W g1t05_2W g1t06_SW
        g1t07_R1 g1t08_R2 g1t09_RS
        g1t10_RW
        g1t11_W1 g1t12_W2 g1t13_WS
        g1t14_WR
        g1t15_S1 g1t16_S2 g1t17_12 g1t18_1S g1t19_21 g1t20_2S
        
        g2t01_1R g2t02_2R g2t03_SR
        g2t04_1W g2t05_2W g2t06_SW
        g2t07_R1 g2t08_R2 g2t09_RS
        g2t10_RW
        g2t11_W1 g2t12_W2 g2t13_WS
        g2t14_WR
        g2t15_S1 g2t16_S2 g2t17_12 g2t18_1S g2t19_21 g2t20_2S
        
        g3t01_1R g3t02_2R g3t03_SR
        g3t04_1W g3t05_2W g3t06_SW
        g3t07_R1 g3t08_R2 g3t09_RS
        g3t10_RW
        g3t11_W1 g3t12_W2 g3t13_WS
        g3t14_WR
        g3t15_S1 g3t16_S2 g3t17_12 g3t18_1S g3t19_21 g3t20_2S

        age sex race smokstatus / covb ties =efron RL=WALD;

    strata type;
    id pptid;
RUN;





* 3 STATE /  6 Transition Types:  multistate survival;
PROC IMPORT OUT= WORK.SAfull
    DATAFILE= "SAfull_allgroups_agesexracesmoking_0306.csv"
    DBMS=CSV REPLACE;
    GETNAMES=YES;
    DATAROW=2;
RUN;



title "SA 3/6";
PROC phreg data = SAfull covs(AGGREGATE) 
    class race(REF=FIRST) smokstatus;
    model time*status(0) = 
        g1t01_NR
        g1t02_NW
        g1t03_RN
        g1t04_RW
        g1t05_WN
        g1t06_WR
        g2t01_NR
        g2t02_NW
        g2t03_RN
        g2t04_RW
        g2t05_WN
        g2t06_WR
        g3t01_NR
        g3t02_NW
        g3t03_RN
        g3t04_RW
        g3t05_WN
        g3t06_WR
        
        age sex race smokstatus / covb ties =efron RL=WALD;

    strata type;
    id pptid;

RUN;



*Bonus:  time varying SA 5/20;
*takes a while to run:  13 hours on 64 bit machine with 16GB RAM;
PROC IMPORT OUT= WORK.SAfull
    DATAFILE= "SAfull_allgroups_agesexracesmoking.csv"
    DBMS=CSV REPLACE;
    GETNAMES=YES;
    DATAROW=2;
RUN;



title "TimeVarying! NO weights;  with covariates";
PROC phreg data = SAfull covs(AGGREGATE) /*outest = cookieEST covout /*noprint*/;
    class race(REF=FIRST) smokstatus;
    model time*status(0) = /*group*/
        g1t01_1R g1t02_2R g1t03_SR
        g1t04_1W g1t05_2W g1t06_SW
        g1t07_R1 g1t08_R2 g1t09_RS
        g1t10_RW
        g1t11_W1 g1t12_W2 g1t13_WS
        g1t14_WR
        g1t15_S1 g1t16_S2 g1t17_12 g1t18_1S g1t19_21 g1t20_2S
        
        g2t01_1R g2t02_2R g2t03_SR
        g2t04_1W g2t05_2W g2t06_SW
        g2t07_R1 g2t08_R2 g2t09_RS
        g2t10_RW
        g2t11_W1 g2t12_W2 g2t13_WS
        g2t14_WR
        g2t15_S1 g2t16_S2 g2t17_12 g2t18_1S g2t19_21 g2t20_2S
        
        g3t01_1R g3t02_2R g3t03_SR
        g3t04_1W g3t05_2W g3t06_`SW
        g3t07_R1 g3t08_R2 g3t09_RS
        g3t10_RW
        g3t11_W1 g3t12_W2 g3t13_WS
        g3t14_WR
        g3t15_S1 g3t16_S2 g3t17_12 g3t18_1S g3t19_21 g3t20_2S
        
        age sex race smokstatus 
        
        timeg1t01_1R 
        timeg1t02_2R 
        timeg1t03_SR 
        timeg1t04_1W 
        timeg1t05_2W 
        timeg1t06_SW 
        timeg1t07_R1 
        timeg1t08_R2 
        timeg1t09_RS 
        timeg1t10_RW 
        timeg1t11_W1 
        timeg1t12_W2 
        timeg1t13_WS 
        timeg1t14_WR 
        timeg1t15_S1 
        timeg1t16_S2 
        timeg1t17_12 
        timeg1t18_1S 
        timeg1t19_21 
        timeg1t20_2S 
        
        timeg2t01_1R 
        timeg2t02_2R 
        timeg2t03_SR 
        timeg2t04_1W 
        timeg2t05_2W 
        timeg2t06_SW 
        timeg2t07_R1 
        timeg2t08_R2 
        timeg2t09_RS 
        timeg2t10_RW 
        timeg2t11_W1 
        timeg2t12_W2 
        timeg2t13_WS 
        timeg2t14_WR 
        timeg2t15_S1 
        timeg2t16_S2 
        timeg2t17_12 
        timeg2t18_1S 
        timeg2t19_21 
        timeg2t20_2S 
        
        timeg3t01_1R 
        timeg3t02_2R 
        timeg3t03_SR 
        timeg3t04_1W 
        timeg3t05_2W 
        timeg3t06_SW 
        timeg3t07_R1 
        timeg3t08_R2 
        timeg3t09_RS 
        timeg3t10_RW 
        timeg3t11_W1 
        timeg3t12_W2 
        timeg3t13_WS 
        timeg3t14_WR 
        timeg3t15_S1 
        timeg3t16_S2 
        timeg3t17_12 
        timeg3t18_1S 
        timeg3t19_21 
        timeg3t20_2S / covb ties =efron RL=WALD;
    strata type;
    id pptid;
    
    timeg1t01_1R = g1t01_1R*log(time); 
    timeg1t02_2R = g1t02_2R*log(time); 
    timeg1t03_SR = g1t03_SR*log(time); 
    timeg1t04_1W = g1t04_1W*log(time); 
    timeg1t05_2W = g1t05_2W*log(time); 
    timeg1t06_SW = g1t06_SW*log(time); 
    timeg1t07_R1 = g1t07_R1*log(time); 
    timeg1t08_R2 = g1t08_R2*log(time); 
    timeg1t09_RS = g1t09_RS*log(time); 
    timeg1t10_RW = g1t10_RW*log(time); 
    timeg1t11_W1 = g1t11_W1*log(time); 
    timeg1t12_W2 = g1t12_W2*log(time); 
    timeg1t13_WS = g1t13_WS*log(time); 
    timeg1t14_WR = g1t14_WR*log(time); 
    timeg1t15_S1 = g1t15_S1*log(time); 
    timeg1t16_S2 = g1t16_S2*log(time); 
    timeg1t17_12 = g1t17_12*log(time); 
    timeg1t18_1S = g1t18_1S*log(time); 
    timeg1t19_21 = g1t19_21*log(time); 
    timeg1t20_2S = g1t20_2S*log(time);
    
    timeg2t01_1R = g2t01_1R*log(time); 
    timeg2t02_2R = g2t02_2R*log(time); 
    timeg2t03_SR = g2t03_SR*log(time); 
    timeg2t04_1W = g2t04_1W*log(time); 
    timeg2t05_2W = g2t05_2W*log(time); 
    timeg2t06_SW = g2t06_SW*log(time); 
    timeg2t07_R1 = g2t07_R1*log(time); 
    timeg2t08_R2 = g2t08_R2*log(time); 
    timeg2t09_RS = g2t09_RS*log(time); 
    timeg2t10_RW = g2t10_RW*log(time); 
    timeg2t11_W1 = g2t11_W1*log(time); 
    timeg2t12_W2 = g2t12_W2*log(time); 
    timeg2t13_WS = g2t13_WS*log(time); 
    timeg2t14_WR = g2t14_WR*log(time); 
    timeg2t15_S1 = g2t15_S1*log(time); 
    timeg2t16_S2 = g2t16_S2*log(time); 
    timeg2t17_12 = g2t17_12*log(time); 
    timeg2t18_1S = g2t18_1S*log(time); 
    timeg2t19_21 = g2t19_21*log(time); 
    timeg2t20_2S = g2t20_2S*log(time);
    
    timeg3t01_1R = g3t01_1R*log(time); 
    timeg3t02_2R = g3t02_2R*log(time); 
    timeg3t03_SR = g3t03_SR*log(time); 
    timeg3t04_1W = g3t04_1W*log(time); 
    timeg3t05_2W = g3t05_2W*log(time); 
    timeg3t06_SW = g3t06_SW*log(time); 
    timeg3t07_R1 = g3t07_R1*log(time); 
    timeg3t08_R2 = g3t08_R2*log(time); 
    timeg3t09_RS = g3t09_RS*log(time); 
    timeg3t10_RW = g3t10_RW*log(time); 
    timeg3t11_W1 = g3t11_W1*log(time); 
    timeg3t12_W2 = g3t12_W2*log(time); 
    timeg3t13_WS = g3t13_WS*log(time); 
    timeg3t14_WR = g3t14_WR*log(time); 
    timeg3t15_S1 = g3t15_S1*log(time); 
    timeg3t16_S2 = g3t16_S2*log(time); 
    timeg3t17_12 = g3t17_12*log(time); 
    timeg3t18_1S = g3t18_1S*log(time); 
    timeg3t19_21 = g3t19_21*log(time); 
    timeg3t20_2S = g3t20_2S*log(time); 
    
*http://www.ats.ucla.edu/stat/examples/asa/test_proportionality.htm;
    proportionality_test: test 
        timeg1t01_1R,
        timeg1t02_2R,
        timeg1t03_SR,
        timeg1t04_1W,
        timeg1t05_2W,
        timeg1t06_SW,
        timeg1t07_R1,
        timeg1t08_R2,
        timeg1t09_RS,
        timeg1t10_RW,
        timeg1t11_W1,
        timeg1t12_W2,
        timeg1t13_WS,
        timeg1t14_WR,
        timeg1t15_S1,
        timeg1t16_S2,
        timeg1t17_12,
        timeg1t18_1S,
        timeg1t19_21,
        timeg1t20_2S,
        timeg2t01_1R,
        timeg2t02_2R,
        timeg2t03_SR,
        timeg2t04_1W,
        timeg2t05_2W,
        timeg2t06_SW,
        timeg2t07_R1,
        timeg2t08_R2,
        timeg2t09_RS,
        timeg2t10_RW,
        timeg2t11_W1,
        timeg2t12_W2,
        timeg2t13_WS,
        timeg2t14_WR,
        timeg2t15_S1,
        timeg2t16_S2,
        timeg2t17_12,
        timeg2t18_1S,
        timeg2t19_21,
        timeg2t20_2S,
        timeg3t01_1R,
        timeg3t02_2R,
        timeg3t03_SR,
        timeg3t04_1W,
        timeg3t05_2W,
        timeg3t06_SW,
        timeg3t07_R1,
        timeg3t08_R2,
        timeg3t09_RS,
        timeg3t10_RW,
        timeg3t11_W1,
        timeg3t12_W2,
        timeg3t13_WS,
        timeg3t14_WR,
        timeg3t15_S1,
        timeg3t16_S2,
        timeg3t17_12,
        timeg3t18_1S,
        timeg3t19_21,
        timeg3t20_2S;
RUN;

