function [ErrorCode] = AnalyseComputationalEnvironment(DataInput, DataInputAdd)
% PARALLEL CONTEXT
% In a parallel context, this function is used to check the cluster defined by the user.
% If no error happen the function returns 0. The function complies with
% Windows/Linux operating systems and Matlab/Octave software.
%
%
% INPUT/OUTPUT description:
%
%
% DataInput
%   is the strcture option_.parallel, with the follow fields:
%
%             Local         1 is on local machine, 0 remote
%      ComputerName         the computer name.
%            CPUnbr         the CPU's
%          UserName         the user name for the ComputerName.
%          Password         the password for the user name in ComputerName.
%       RemoteDrive         Drive used for Remote computation (data exchange, etc): must be contain 'RemoteFolder'.
%   RemoteDirectory         Folder in RemoteDrive used for Remote computation.
%  MatlabOctavePath         Path to MATLAB or Octave executable.
%        DynarePath         Path to matlab directory within the Dynare installation directory.
%
%   This information is typed by the user in the DYNARE configuration file and is parsed by the preprocessor,
%   the goal of this function is to check if configuration is correct and if dynare
%   can be executed successfully in parallel mode.
%
%
% DataInputAdd
%   it is the structure options_.parallel_info. Currently , only the string in the
%   field RemoteTmpFolder (the temporary directory created/destroyed on remote
%   computer) is used.

if ispc, 
    [tempo, MasterName]=system('hostname');
    MasterName=deblank(MasterName);
end

RemoteTmpFolder=DataInputAdd.RemoteTmpFolder;
dynareParallelMkDir(RemoteTmpFolder,DataInput);


% The variable ErrorCode is initialized at 0. If there are non problems with
% Local, ComputerName connections,... in general with parallel software execution,
% the ErrorCode is unchanged, in the others cases 1, 2 , ... The values
% table is below.
%
%
%   Table for ErrorCode Values.
%
%   ErrorCode -> 0      Initial Value -> No Error Detected!!!
%   ErrorCode -> 1 ...  When an error is detected, the values 1, 2, 3... are
%   used to specify the type of error or warning.
%
%   Value 1:    The variable 'Local' has a bad value!
%
%   Value 2:    The variable 'CPUnbr' has a bad value. For more information
%               see http://www.dynare.org/DynareWiki/ParallelDynare.
%         2.1   [warning] The user asks to use more CPU's than those available.
%         2.2   [warning] There are unused CPU's!
%
%
%   Value 3:    The remote computer is unreachable!!!
%
%   Value 4:    The fields user name and/or password are/is empty!
%
%   Value 5:    Remote Drive and/or Remote Folder do not exist!
%
%   Value 6:    It is impossible write/read files on the remote computer.
%
%   Value 7:    The values user and/or passwd are incorrect or the user has
%               no permissions to execute a Matlab session. Or simply
%               Matlab path (MatlabOctavePath) is incorrect!
%
%   Value 8:    Dynare path (DynarePath) is incorrect!
%
%   Value 9:    It is impossible delete remote computational temporary files!
%
%
%
%
% Currently when errors are detected execution simply stops and users can
% fix configuration errors according to the error type.

% Copyright (C) 2009-2012 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.


ErrorCode=0;


for Node=1:length(DataInput) % To obtain a recoursive function remove the 'for'
    % and use AnalyseComputationalEnvironment with differents input!
    
    
    % Determine the operating system or software version when necessary
    % for different command types.
   
    OScallerUnix=~ispc;
    OScallerWindows=ispc;
    OStargetUnix=strcmpi('unix',DataInput(Node).OperatingSystem);
    if isempty(DataInput(Node).OperatingSystem),
        OStargetUnix=OScallerUnix;
    end
    OStargetWindows=strcmpi('windows',DataInput(Node).OperatingSystem);
    if isempty(DataInput(Node).OperatingSystem),
        OStargetWindows=OScallerWindows;
    end
    
    Environment= (OScallerUnix || OStargetUnix);
    
    disp(' ');
    disp(' ');
    disp(['Testing computer -> ',DataInput(Node).ComputerName,' <- ...']);
    disp(' ');
    disp(' ');
    
    % The function is composed by two main blocks, determined by the 'Local'
    % variable.
    
    % This check can be removed ... according to the dynare parser
    % strategy.
    
    if ((DataInput(Node).Local == 0) |(DataInput(Node).Local == 1))
        % Continue it is Ok!
        disp('Check on Local Variable ..... Ok!');
        disp(' ');
        disp(' ');
        
    else
        disp('The variable "Local" has a bad value!');
        disp(' ');
        disp('ErrorCode 1.');
        disp(' ');
        disp(' ');
        ErrorCode=1;
        return
        
    end
    
    %         %%%%%%%%%%  Local (No Network) Computing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         Here only the multi-core, or multi-processor avaiable on local
    %         machine are involved in parallel computing. No network
    %         comunications are required!
    
    
    % In this case we need to check only the variable 'CPUnbr'.
    
    % We run the parallel code on local computer, so the others fields are automatically
    % fixed by Dynare parser. Then the user can also fill them with wrong values.
    
    
    %         %%%%%%%%%%  Cluster Computing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         Here we can have many computer with multi-core, or multi-processor avaiable on the
    %         network and involved in parallel computing.
    %         So in this case we need more sophisticated check.
    
    
    if (DataInput(Node).Local == 0)
        
        % Now we verify if it is possibile to be connected with the
        % remote computer.
        
        si1=[];
        de1=[];
        
        if Environment
            if OScallerWindows
                [si1 de1]=system(['ping ', DataInput(Node).ComputerName]);
            else
                [si1 de1]=system(['ping ', DataInput(Node).ComputerName, ' -c 4']);
            end
        else
            [si1 de1]=system(['ping ', DataInput(Node).ComputerName]);
        end
        
        if (si1)
            disp(['It is impossibile to be connected to the computer with name "',DataInput(Node).ComputerName,'" using the network!']);
            disp(' ');
            disp('ErrorCode 3.');
            ErrorCode=3;
            disp(' ');
            disp(' ');
            return;
        else
            disp('Check on ComputerName Variable ..... Ok!');
            disp(' ');
            disp(' ');
        end
        
        
        % Now we verify if user name and password are correct and if remote
        % drive and remote folder exist on the remote computer and it is
        % possible to exchange data with them.
        
        if Environment
            % This check can be removed ... according to the dynare parser
            % strategy.
            
            if (isempty(DataInput(Node).UserName))
                disp('The fields UserName is empty!');
                disp(' ');
                disp('ErrorCode 4.');
                disp(' ');
                disp(' ');
                ErrorCode=4;
                return
            end
            disp('Check on UserName Variable ..... Ok!');
            disp(' ');
            disp(' ');
            
            % This check can be removed ... according to the dynare parser
            % strategy.           
                if (~isempty(DataInput(Node).Password))
                    disp('[WARNING] The field Password should be empty under unix or mac!');
                    disp(' ');
                    disp(['Remove the string ',DataInput(Node).Password,' from this field!']);
                    disp(' ');
                    disp('ErrorCode 4.');
                    disp(' ');
                    disp(' ');
                    ErrorCode=4;
                    %                 return
                else
                    disp('Check on Password Variable ..... Ok!');
                    disp(' ');
                    disp(' ');
                end
        else
            
            % This check can be removed ... according to the dynare parser
            % strategy.
            
            if (isempty(DataInput(Node).UserName)) || (isempty(DataInput(Node).Password))
                disp('The fields UserName and/or Password are/is empty!');
                disp(' ');
                disp('ErrorCode 4.');
                disp(' ');
                disp(' ');
                ErrorCode=4;
                return
            end
            disp('Check on UserName Variable ..... Ok!');
            disp(' ');
            disp(' ');
            disp('Check on Password Variable ..... Ok!');
            disp(' ');
            disp(' ');
            
        end
        
        % Now we very if RemoteDrive and/or RemoteDirectory exist on remote
        % computer!
        
        if Environment
            
            % This check can be removed ... according to the dynare parser
            % strategy.
            
            if  isempty(DataInput(Node).RemoteDirectory)
                disp('The field RemoteDirectory is empty!');
                disp(' ');
                disp('ErrorCode 5.');
                disp(' ');
                disp(' ');
                ErrorCode=5;
                return
            end
            
            % This check can be removed ... according to the dynare parser
            % strategy.
            
                if (~isempty(DataInput(Node).RemoteDrive))
                    disp('[WARNING] The fields RemoteDrive should be empty under unix or mac!');
                    disp(' ');
                    disp(['remove the string ',DataInput(Node).RemoteDrive,' from this field!']);
                    disp(' ');
                    disp('ErrorCode 5.');
                    disp(' ');
                    disp(' ');
                    ErrorCode=5;
                    %                 return
                end
            
            si2=[];
            de2=[];
            if ~isempty(DataInput(Node).Port),
                ssh_token = ['-p ',DataInput(Node).Port];
            else
                ssh_token = '';
            end
            
            [si2 de2]=system(['ssh ',ssh_token,' ',DataInput(Node).UserName,'@',DataInput(Node).ComputerName,' ls ',DataInput(Node).RemoteDirectory,'/',RemoteTmpFolder,'/']);
            
            if (si2)
                disp ('Remote Directory does not exist or is not reachable!');
                disp(' ');
                disp('ErrorCode 5.');
                disp(' ');
                disp(' ');
                ErrorCode=5;
                return
            end
            
            disp('Check on RemoteDirectory Variable ..... Ok!');
            disp(' ');
            disp(' ');
            disp('Check on RemoteDrive Variable ..... Ok!');
            disp(' ');
            disp(' ');
            
        else
            % This check can be removed ... according to the dynare parser
            % strategy.
            
            if (isempty(DataInput(Node).RemoteDrive)||isempty(DataInput(Node).RemoteDirectory))
                disp('Remote RemoteDrive and/or RemoteDirectory is/are empty!');
                disp(' ');
                disp('ErrorCode 5.');
                disp(' ');
                disp(' ');
                ErrorCode=5;
                return
            end
            
            
            si2=[];
            de2=[];
            [si2 de2]=system(['dir \\',DataInput(Node).ComputerName,'\',DataInput(Node).RemoteDrive,'$\',DataInput(Node).RemoteDirectory,'\',RemoteTmpFolder]);
            
            if (si2)
                disp ('Remote Directory does not exist or it is not reachable!');
                disp(' ');
                disp('ErrorCode 5.');
                disp(' ');
                disp(' ');
                ErrorCode=5;
                return
            end
            
            disp('Check on RemoteDirectory Variable ..... Ok!');
            disp(' ');
            disp(' ');
            disp('Check on RemoteDrive Variable ..... Ok!');
            disp(' ');
            disp(' ');
            
        end
        
        
        % Now we verify if it possible to exchange data with the remote
        % computer:
        
        
        % Build a command file to test the matlab execution and dynare path ...
        
        fid = fopen('Tracing.m', 'w+');
        s1=(['fT = fopen(''MatlabOctaveIsOk.txt'',''w+'');\n']);
        s2='fclose(fT);\n';
        SBS=strfind(DataInput(Node).DynarePath,'\');
        DPStr=DataInput(Node).DynarePath;
        if isempty(SBS),
            DPStrNew=DPStr;
        else
            DPStrNew=[DPStr(1:SBS(1)),'\'];
            for j=2:length(SBS),
                DPStrNew=[DPStrNew,DPStr(SBS(j-1)+1:SBS(j)),'\'];
            end
            DPStrNew=[DPStrNew,DPStr(SBS(end)+1:end)];
        end
        s3=['addpath(''',DPStrNew,'''),\n'];
        s4=['try,\n  dynareroot = dynare_config();\n'];
        s41=(['  fT = fopen(''DynareIsOk.txt'',''w+'');\n']);
        s42='  fclose(fT);\n';
        s5=['catch,end,\n'];
        s51=(['  fT = fopen(''DynareFailed.txt'',''w+'');\n']);
        s52='  fclose(fT);\n';
        send='exit';
        StrCommand=([s1,s2,s3,s4,s41,s42,s5,s51,s52,send]);
        
        % Mettere controllo su NbW ...
        % if exist('OCTAVE_VERSION')
        %     NbW = fprintf(fid,StrCommand, '%s');
        % else
        NbW = fprintf(fid,StrCommand, '%s');
        % end
        fclose(fid);
        
        dynareParallelSendFiles('Tracing.m', RemoteTmpFolder,DataInput(Node));
        FindTracing = dynareParallelDir('Tracing.m', RemoteTmpFolder,DataInput(Node));
        
        delete ('Tracing.m');
        
        if (isempty(FindTracing))
            disp ('It is impossible to exchange data with Remote Drive and/or Remote Directory! ErrorCode 6.');
            disp(' ');
            disp('ErrorCode 6.');
            disp(' ');
            disp(' ');
            ErrorCode=6;
            return
        else
            disp('Check on Exchange File with Remote Computer ..... Ok!');
            disp(' ');
            disp(' ');
        end
        
        
        % Now we verify if it is possible execute a matlab/octave section on remote
        % machine when the user is .UserName with password .Password and
        % the path is MatlabOctavePath.
        
        if Environment,
            if ~isempty(DataInput(Node).Port),
                ssh_token = ['-p ',DataInput(Node).Port];
            else
                ssh_token = '';
            end
            if strfind([DataInput(Node).MatlabOctavePath], 'octave') % Hybrid computing Matlab(Master)->Octave(Slaves) and Vice Versa!
                system(['ssh ',ssh_token,' ',DataInput(Node).UserName,'@',DataInput(Node).ComputerName,' "cd ',DataInput(Node).RemoteDirectory,'/',RemoteTmpFolder,  '; ', DataInput(Node).MatlabOctavePath, ' Tracing.m;" &']);
            else
                system(['ssh ',ssh_token,' ',DataInput(Node).UserName,'@',DataInput(Node).ComputerName,' "cd ',DataInput(Node).RemoteDirectory,'/',RemoteTmpFolder,  '; ', DataInput(Node).MatlabOctavePath, ' -nosplash -nodesktop -minimize -r Tracing;" &']);
            end
        else
            if ~strcmp(DataInput(Node).ComputerName,MasterName), % run on remote machine
                if  strfind([DataInput(Node).MatlabOctavePath], 'octave') % Hybrid computing Matlab(Master)->Octave(Slaves) and Vice Versa!
                    [NonServeS NenServeD]=system(['start /B psexec \\',DataInput(Node).ComputerName,' -e -u ',DataInput(Node).UserName,' -p ',DataInput(Node).Password,' -W ',DataInput(Node).RemoteDrive,':\',DataInput(Node).RemoteDirectory,'\',RemoteTmpFolder ' -low   ',DataInput(Node).MatlabOctavePath,' Tracing.m']);
                else
                    [NonServeS NenServeD]=system(['start /B psexec \\',DataInput(Node).ComputerName,' -e -u ',DataInput(Node).UserName,' -p ',DataInput(Node).Password,' -W ',DataInput(Node).RemoteDrive,':\',DataInput(Node).RemoteDirectory,'\',RemoteTmpFolder ' -low   ',DataInput(Node).MatlabOctavePath,' -nosplash -nodesktop -minimize -r Tracing']);
                end
            else % run on local machine via the network: user and passwd cannot be used!
                if  strfind([DataInput(Node).MatlabOctavePath], 'octave') % Hybrid computing Matlab(Master)->Octave(Slaves) and Vice Versa!
                    [NonServeS NenServeD]=system(['start /B psexec \\',DataInput(Node).ComputerName,' -e ',' -W ',DataInput(Node).RemoteDrive,':\',DataInput(Node).RemoteDirectory,'\',RemoteTmpFolder ' -low   ',DataInput(Node).MatlabOctavePath,' Tracing.m']);
                else
                    [NonServeS NenServeD]=system(['start /B psexec \\',DataInput(Node).ComputerName,' -e ',' -W ',DataInput(Node).RemoteDrive,':\',DataInput(Node).RemoteDirectory,'\',RemoteTmpFolder ' -low   ',DataInput(Node).MatlabOctavePath,' -nosplash -nodesktop -minimize -r Tracing']);
                end
            end
                
        end
        
        % Timer da fissare, nei valori di attesa!
        
        t1=fix(clock);
        
        if t1(5)+1>60;
            t2=2;
        else t2=t1(5)+1;
        end
        
        Flag=0;
        
        while (1);
            if Flag==0
                disp('Try to run matlab/octave on remote machine ... ');
                disp(' ');
                disp('please wait ... ');
                disp(' ');
                Flag=1;
            end
            nt=fix(clock);
            nt(5)-t2;
            
            if (~isempty (dynareParallelDir('MatlabOctaveIsOk.txt',RemoteTmpFolder,DataInput(Node)))) || ((nt(5)-t2)>0)
                if ((nt(5)-t2)>0)
                    ErrorCode=7;
                end
                break
            end
            
        end
        
        if  (ErrorCode==7)
            
            disp ('It is not possible execute a matlab session on remote machine!');
            disp(' ');
            disp('ErrorCode 7.');
            disp(' ');
            disp(' ');
            ErrorCode=7;
            dynareParallelRmDir(RemoteTmpFolder,DataInput(Node));
            return
            
        else
            disp('Check on MatlabOctave Path and MatlabOctave Program Execution on remote machine ..... Ok!');
            disp(' ');
            disp(' ');
            
            % Now we verify if the DynarePath is correct ...
            disp('Check the Dynare path on remote machine ... ');
            disp(' ');
            disp('please wait ... ');
            disp(' ');
            pause(2)
            
            if isempty(dynareParallelDir('DynareIsOk.txt',RemoteTmpFolder,DataInput(Node)))
                ErrorCode=8;
            end
            
            if  (ErrorCode==8)
                
                disp ('The DynarePath is incorrect!');
                disp(' ');
                disp('ErrorCode 8.');
                disp(' ');
                disp(' ');
                ErrorCode=8;
                dynareParallelRmDir(RemoteTmpFolder,DataInput(Node));
                return
                
            else
                disp('Check on Dynare Path remote machine ..... Ok!');
                disp(' ');
                disp(' ');
            end
        end
        
        
        % Now we verify if it is possible delete remote computational traces!
        
        dynareParallelRmDir(RemoteTmpFolder,DataInput(Node));
        
        si3=[];
        
        si3=dynareParallelDir('Tracing.m', RemoteTmpFolder,DataInput(Node));
        
        if (isempty(si3))
            disp ('Check on Delete Remote Computational Traces ..... Ok!');
            disp(' ');
            disp(' ');
        else
            disp ('It is impossible to delete temporary files on remote machine!');
            disp(' ');
            disp('ErrorCode 9.');
            disp(' ');
            disp(' ');
            ErrorCode=9;
            return
        end
        
        
    end
    % Now we check the variable 'CPUnbr'.
    
    % This check can be removed ... according to the dynare parser
    % strategy.
    
    yn=isempty(DataInput(Node).CPUnbr);
    
    if yn==1
        % The field is empty!
        disp('The field "CPUnbr" is empty!');
        disp(' ');
        disp('ErrorCode 2.');
        disp(' ');
        disp(' ');
        ErrorCode=2;
        return
    end
    
    % This check can be removed ... according to the dynare parser
    % strategy.
    
    
    
    % We look for the information on local computer hardware.
    
    si0=[];
    de0=[];
    
    disp('Checking Hardware please wait ...');
    if (DataInput(Node).Local == 1)
        if Environment,
            [si0 de0]=system('grep processor /proc/cpuinfo');
        else
            [si0 de0]=system(['psinfo \\']);
        end
    else
        if Environment,
            if ~isempty(DataInput(Node).Port),
                ssh_token = ['-p ',DataInput(Node).Port];
            else
                ssh_token = '';
            end
            if OStargetUnix,
                [si0 de0]=system(['ssh ',ssh_token,' ',DataInput(Node).UserName,'@',DataInput(Node).ComputerName,' grep processor /proc/cpuinfo']);
            else
                [si0 de0]=system(['ssh ',ssh_token,' ',DataInput(Node).UserName,'@',DataInput(Node).ComputerName,' psinfo']);
            end
        else
            [si0 de0]=system(['psinfo \\',DataInput(Node).ComputerName,' -u ',DataInput(Node).UserName,' -p ',DataInput(Node).Password]);
        end
    end
    
    
    RealCPUnbr='';
%    keyboard;
    RealCPUnbr=GiveCPUnumber(de0);
    
    % Questo controllo penso che si possa MIGLIORARE!!!!!
    if  isempty (RealCPUnbr) && Environment==0,
        [si0 de0]=system(['psinfo \\',DataInput(Node).ComputerName]);
    end        
    RealCPUnbr=GiveCPUnumber(de0);

    if  isempty (RealCPUnbr)
        % An error occurred when we try to know the Cpu/Cores
        % numbers.
        disp('It is impossible determine the number of Cpu/Processor avaiable on this machine!');
        disp(' ');
        disp('ErrorCode 2.');
        disp(' ');
        if Environment
            disp('Check the command "$less /proc/cpuinfo" ... !');
        else
            disp('Check if the pstools are installed and are in machine path! And check the command "psinfo \\"');
        end
        disp(' ');
        ErrorCode=2;
        return
    end
    
    
    % Trasforming the input data provided in a form [n1:n2] in a single numerical
    % value.
    
    
    CPUnbrUser=length(DataInput(Node).CPUnbr);
    maxCPUnbrUser=max(DataInput(Node).CPUnbr)+1;
    
    disp(['Hardware has ', num2str(RealCPUnbr),' Cpu/Cores!']);
    disp(['User requires ',num2str(CPUnbrUser),' Cpu/Cores!']);
    if  CPUnbrUser==RealCPUnbr,
        % It is Ok!
        disp('Check on CPUnbr Variable ..... Ok!');
        disp(' ');
        disp(' ');
        disp(' ');
        
    end
    
    if CPUnbrUser > RealCPUnbr
        disp('Warning! The user asks to use more CPU''s than those available.');
        disp(' ');
        disp(' ');
        ErrorCode=2.1;
        % return
        
    end
    if CPUnbrUser < RealCPUnbr
        disp('Warning! There are unused CPU''s!');
        disp(' ');
        disp(' ');
        ErrorCode=2.2;
        % return
    end
    
    disp(['Test for Cluster computation, computer ',DataInput(Node).ComputerName, ' ..... Passed!']);
    disp(' ');
    disp(' ');
    
    
end

return