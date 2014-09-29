
for test =1:5
    
    RootName = sprintf ('Test%02.0f', test); 
    FileName = [ RootName '.SD.out' ];
    PlotFASToutput( {[RootName filesep 'NREL_Results' filesep FileName],[RootName filesep FileName]}, ...
                    {'NREL results','CertTest'});
                
%                 'M3J2MKze'
%                 'M3J1MKze'
%                 'M3J1MKye'
end