% Select a folder with calibrated and modulated spectra
	%       Version 4 is meant to reduce memory usage by deleting SpectraDataMod info and working toward slope
	%       calculation and diffusivity calculation
	%       Version 5 includes NN application for temperature
	%%
	%   1 - The first part finds all the calibration files which will then each be
	%       processed and joined together to make a theta to train the NN with.
	%                   DONE
	%   2 - Find the green, fluoro, and bump indices based on the lowest
	%       temperature and closest to the origin file
    %   3 - The NN is trained to based on PeakIntensity, IntegratedIntensity,
	%       FWHM, PeakWL, Ratio, and Normalized Locations from theta.  
	%   4 - Then the modulated spectra at each temperature will be processed with a
	%       theta created for each motor position.
	%   5 - The NN is applied to each motor position to get a temperature.
	%   6 - The fft of the reference signal is then determined, followed by the fft
	%       of each spectra characteristic, followed by the NN temperature, 
	%       and then the ln(MAG(T)) and PHASE(T) at the modulation frequency at
	%       each position is stored for Peak, Integrated, PeakWL, and NN temperature.
	%   7 - The slope of an absolute value function is then fit to the data for the
	%       MAG and PHASE.
	%   8 - The diffusivity is then calculated as the [-Pi*Freq]/[slope(MAG)xslope(PHASE)]
	%   
	clc
	% close all
26	clear all
27	
28	%%
29	[ Folder ] = SharedSpectraFolders;
30	cd(Folder.Input)
31	open('LeuvenSpectralSupportProgramInputParameters.m')
32	pause
33	[ Inputs ] = LeuvenSpectralSupportProgramInputParameters;
34	
35	Folder.SpectraExpHome=uigetdir([Folder.Doc,Folder.mainDir,Folder.fspec]);
36	cd(Folder.SpectraExpHome)
37	dirs = regexp(genpath(pwd),['[^;]*'],'match');
38	for i = 2:length(dirs)
39	    Folder.TempCell{1,i-1}=dirs{1,i};
40	end
41	xstruct=dir;
42	xcell=struct2cell(xstruct);
43	for i = 1:length(xstruct)
44	    xcellnames{1,i}=xcell{1,i};
45	end
46	temp = regexp(xcellnames,'.*\.calib','match');
47	temp(cellfun(@isempty,temp)) = [];
48	for i = 1:length(temp)
49	    FileName.CalibCell{1,i} = temp{1,i}{1,1};
50	end
51	
52	%% Try and find the Indices for Fluorescence, Grn, Bump, etc.
53	[y,Fs] = audioread([Folder.Sound,'\Finding Indices.mp3']);sound(y,Fs);
54	cd(Folder.Programs)
55	LeuvenSpectralAnalysisProgramFindFluoroIndices_V2;
56	
57	%% 1 - Calibration files processing
58	FileName.numCalibCell = length(FileName.CalibCell);
59	FileName.SpecExpHomeFolder = Folder.SpectraExpHome;
60	tempAll_1  = [];
61	tempNorm_1 = [];
62	[y,Fs] = audioread([Folder.Sound,'\Calibrating Files.mp3']);sound(y,Fs);
63	
64	for CalibI = 1:FileName.numCalibCell    
65	    FileName.TxtName = FileName.CalibCell{1,CalibI};
66	    %Save current FileName without .calib at end
67	    FileName.FullName = FileName.TxtName(1:end-6);
68	    %Save full path of current FileName 
69	    FileName.FullPath = horzcat(FileName.SpecExpHomeFolder,'\',FileName.TxtName);
70	    %Save Time Stamp of current FileName
71	    scrap = dir(FileName.FullPath);
72	    FileName.Time=scrap.date;    
73	    FileName.TimeNum=scrap.datenum;
74	%     FileName.DateStr=FileName.Time(1:6);
75	    FileName.DateStr=[FileName.DateExp(4:5),'-',FileName.DateExp(1:3)];
76	    FileName.YearStr=scrap.date(8:11);    clear scrap
77	    cd(Folder.Programs)
78	    [SpectraData] = LeuvenSpectralAnalysisProgramSpecCalib_V2(Folder,FileName,Inputs,IndicesFluoro,IndicesGrn,IndicesBump,FileData);
79	    
80	    tempAll_2  = SpectraData.ThetaAll;
81	    tempAll_3  = [tempAll_1;tempAll_2];
82	    tempAll_1  =  tempAll_3;
83	        
84	    tempNorm_2  = SpectraData.ThetaNorm;
85	    tempNorm_3  = [tempNorm_1;tempNorm_2];
86	    tempNorm_1  =  tempNorm_3;
87	    
88	    TempTC(:,CalibI) = SpectraData.TC;
89	    TempK(:,CalibI)  = SpectraData.TempK;
90	    TempC(:,CalibI)  = SpectraData.TempC;
91	    
92	end
93	SpectraData.ThetaAll  = tempAll_3;
94	SpectraData.ThetaNorm = tempNorm_3;
95	
96	cd(Folder.Temperature)
97	TempFolderName = [FileData.monthStr,FileData.dateStr];
98	mkdir(TempFolderName);
99	cd([Folder.Temperature,'\',TempFolderName])
100	dlmwrite([TempFolderName,', TC.txt'], TempTC)
101	dlmwrite([TempFolderName,', TempK.txt'], TempK)
102	dlmwrite([TempFolderName,', TempC.txt'], TempC)
103	
104	%% Run NN Calibration Process
105	[y,Fs] = audioread([Folder.Sound,'\Training Neural Network.mp3']);sound(y,Fs);
106	cd(Folder.Programs)
107	LeuvenSpectralAnalysisProgramNNCalib_V1;
108	
109	%% Process Modulated Spectra
110	%remove certain SpectraDataMod fields to reduce memory usage
111	fields = {'BaseValues','Raw','BaseSubtracted','SmoothedNorm','Bump','SumBump',...
112	    'TC','TempC','TempK','Smoothed','ThetaNorm'};
113	[y,Fs] = audioread([Folder.Sound,'\Smoothing Spectra.mp3']);
114	[y2,Fs2] = audioread([Folder.Sound,'\Done Processing Spectra.mp3']);
115	
116	for FolderLoop = 1:length(Folder.TempCell)
117	    sound(y,Fs);
118	    %Counts the number of 
119	    files = dir(Folder.TempCell{1,FolderLoop}); % equal
120	    files(strncmp({files.name}, '.', 1)) = []; % new, no exceptions
121	    
122	    for FileLoop = 1:length(files)
123	        cd(Folder.Programs)
124	        FileName.FullPath = [Folder.TempCell{1,FolderLoop},'\',files(FileLoop,1).name];
125	        FileName.CurrentFileFullPath = FileName.FullPath;
126	        
127	        [extra] = LeuvenSpectralAnalysisProgramSpecModulated_V1...
128	            (Folder,FileName,Inputs,IndicesFluoro,IndicesGrn,IndicesBump,FolderLoop);
129	%         extra = rmfield(extra,fields);
130	        SpectraDataMod(FolderLoop,FileLoop) = extra;
131	%         Spectra temp too
132	        clear extra
133	%         SpectraDataMod = rmfield(SpectraDataMod,fields);
134	        close all
135	    end
136	end
137	sound(y2,Fs2);
138	fprintf('            Done Smoothing Spectra')
139	
140	%%  Apply NN to each motorposition
141	%       remove certain SpectraDataMod fields to reduce memory usage
142	% fields = {'ThetaNorm'};   %'ThetaAll',
143	clear extra
144	for FolderLoop = 1:length(Folder.TempCell)
145	    files = dir(Folder.TempCell{1,FolderLoop}); % equal
146	    files(strncmp({files.name}, '.', 1)) = []; % new, no exceptions
147	    for FileLoop = 1:length(files)
148	        SpectraDataMod(FolderLoop,FileLoop).NNTemp={[],[],[],[]};
149	        SpectraDataMod(FolderLoop,FileLoop).ThetaNorm=[];
150	    end
151	end
152	[y,Fs] = audioread([Folder.Sound,'\Applying Neural Network.mp3']);sound(y,Fs);
153	for FolderLoop = 1:length(Folder.TempCell)
154	    %Counts the number of 
155	    files = dir(Folder.TempCell{1,FolderLoop}); % equal
156	    files(strncmp({files.name}, '.', 1)) = []; % new, no exceptions
157	    
158	    
159	    for FileLoop = 1:length(files)
160	                
161	        cd(Folder.Programs)
162	        FileName.FullPath = [Folder.TempCell{1,FolderLoop},'\',files(FileLoop,1).name];
163	        FileName.CurrentFileFullPath = FileName.FullPath;
164	        
165	        [extra, NNs] = LeuvenSpectralAnalysisProgramNNApplication_V1...
166	            (Folder,FileName,IndicesFluoro,SpectraDataMod,FolderLoop,FileLoop);
167	        SpectraDataMod(FolderLoop,FileLoop) = extra;
168	        
169	%         extra = rmfield(extra,fields);
170	        clear extra
171	    end
172	end
173	
174	%% Perform FFT
175	[y,Fs] = audioread([Folder.Sound,'\Perform FFT.mp3']);sound(y,Fs);
176	% Folder loop is temperature indices
177	for FolderLoop = 1:length(Folder.TempCell)
178	    %Counts the number of 
179	    files = dir(Folder.TempCell{1,FolderLoop}); % equal
180	    files(strncmp({files.name}, '.', 1)) = []; % new, no exceptions
181	    
182	    for FileLoop = 1:length(files)
183	        
184	        %Temperature storage
185	        SpectraDataMod(FolderLoop,FileLoop).Temperature = Folder.TempCell{FolderLoop}(end-3:end-1);
186	        SpectraDataMod(FolderLoop,FileLoop).TempStr = Folder.TempCell{FolderLoop}(end-3:end-1);
187	        SpectraDataMod(FolderLoop,FileLoop).TempNum = str2double(SpectraDataMod(FolderLoop,FileLoop).TempStr);
188	        
189	        %Perform FFT
190	        SpectraDataMod(FolderLoop,FileLoop).fftRefGrn      = ...
191	            fft(SpectraDataMod(FolderLoop,FileLoop).SumGrn);
192	        SpectraDataMod(FolderLoop,FileLoop).fftPeak        = ...
193	            fft(SpectraDataMod(FolderLoop,FileLoop).PeakIntensity);
194	        SpectraDataMod(FolderLoop,FileLoop).fftIntegrated   = ...
195	            fft(SpectraDataMod(FolderLoop,FileLoop).IntegIntensity);
196	        SpectraDataMod(FolderLoop,FileLoop).fftPeakWL      = ...
197	            fft(SpectraDataMod(FolderLoop,FileLoop).PeakWL);
198	        for iNNs = 1:NNs
199	            SpectraDataMod(FolderLoop,FileLoop).fftNNTemp{1,iNNs}  = ...
200	            fft(SpectraDataMod(FolderLoop,FileLoop).NNTemp{1,iNNs});
201	        end
202	        
203	                %% spacing of t
204	                t=SpectraDataMod(FolderLoop,FileLoop).Time;
205	                for i=2:length(t)
206	                    dt(i-1)=t(i)-t(i-1);
207	                end
208	                N=length(SpectraDataMod(FolderLoop,FileLoop).Time);
209	                N2=ceil(N/2);
210	        SpectraDataMod(FolderLoop,FileLoop).dt=mean(dt);
211	                Fs=1/SpectraDataMod(FolderLoop,FileLoop).dt;
212	        SpectraDataMod(FolderLoop,FileLoop).fftFrequencies = 0:Fs/length(t):Fs/2;
213	        
214	        %Find modulation frequency based on SpectraDataMod.fftRefGrn
215	            [pks,locs] = findpeaks((abs(SpectraDataMod(FolderLoop,FileLoop).fftRefGrn)),...
216	                'MINPEAKDISTANCE',10);%,'THRESHOLD',1,'NPEAKS',2)
217	            [t,b] = max(pks);
218	            ModFreqIndex = locs(b);
219	            ModFreq = SpectraDataMod(FolderLoop,FileLoop).fftFrequencies(ModFreqIndex);
220	        SpectraDataMod(FolderLoop,FileLoop).ModFreqIndex = ModFreqIndex;
221	        SpectraDataMod(FolderLoop,FileLoop).ModFreq = ModFreq;
222	        
223	        
224	%         %Plot Phase of reference signal, Peak, Integ, PeakWL, and NN
225	%         figure;  hold all
226	%         plot(SpectraDataMod(FolderLoop,FileLoop).fftFrequencies,...
227	%             unwrap(angle(SpectraDataMod(FolderLoop,FileLoop).fftRefGrn(1:N2,1))));
228	%         plot(SpectraDataMod(FolderLoop,FileLoop).fftFrequencies,...
229	%             unwrap(angle(SpectraDataMod(FolderLoop,FileLoop).fftPeak(1:N2,1))));
230	%         plot(SpectraDataMod(FolderLoop,FileLoop).fftFrequencies,...
231	%             unwrap(angle(SpectraDataMod(FolderLoop,FileLoop).fftIntegrated(1:N2,1))));
232	%         plot(SpectraDataMod(FolderLoop,FileLoop).fftFrequencies,...
233	%             unwrap(angle(SpectraDataMod(FolderLoop,FileLoop).fftPeakWL(1:N2,1))));
234	% %         plot(SpectraDataMod(FolderLoop,FileLoop).fftFrequencies,...
235	% %             unwrap(angle(SpectraDataMod(FolderLoop,FileLoop).fftNNTemp(1:N2,1))));
236	%         legend('Reference','Peak','Integ','Peak WL','NN Temp')
237	%         Pos     = num2str(SpectraDataMod(FolderLoop,FileLoop).MotorPosition);
238	%         Temp    = SpectraDataMod(FolderLoop,FileLoop).TempStr;
239	%         TempPos = [Temp,'K, Position ',Pos,' mm'];
240	%         title(TempPos)        
241	%         xlabel('Frequencies (Hz)','fontsize',20)
242	%         ylabel('Phase (rads)','fontsize',20)
243	%         axis('tight') 
244	%         h10=gcf;
245	%         set(gca,'fontsize',15)
246	%         saveas(h10,['(',FileName.DateStr,')Phase, ',TempPos],'jpg')
247	%         hold off
248	%         
249	        
250	        
251	        %Find values at lock-in frequency
252	        SpectraDataMod(FolderLoop,FileLoop).LockInPeak = ...
253	            SpectraDataMod(FolderLoop,FileLoop).fftPeak(ModFreqIndex);
254	        SpectraDataMod(FolderLoop,FileLoop).LockInIntegrated = ...
255	            SpectraDataMod(FolderLoop,FileLoop).fftIntegrated(ModFreqIndex);
256	        SpectraDataMod(FolderLoop,FileLoop).LockInPeakWL = ...
257	            SpectraDataMod(FolderLoop,FileLoop).fftPeakWL(ModFreqIndex);
258	        for iNNs = 1:NNs
259	            SpectraDataMod(FolderLoop,FileLoop).LockInNNTemp{1,iNNs}  = ...
260	            SpectraDataMod(FolderLoop,FileLoop).fftNNTemp{1,iNNs}(ModFreqIndex);
261	        end
262	
263	        %MAGnitude
264	        SpectraDataMod(FolderLoop,FileLoop).LnMagPeak = ...
265	            log(abs(SpectraDataMod(FolderLoop,FileLoop).LockInPeak));
266	        SpectraDataMod(FolderLoop,FileLoop).LnMagIntegrated = ...
267	            log(abs(SpectraDataMod(FolderLoop,FileLoop).LockInIntegrated));
268	        SpectraDataMod(FolderLoop,FileLoop).LnMagPeakWL = ...
269	            log(abs(SpectraDataMod(FolderLoop,FileLoop).LockInPeakWL));
270	        for iNNs = 1:NNs
271	            SpectraDataMod(FolderLoop,FileLoop).LnMagNNTemp{1,iNNs}  = ...
272	            log(abs(SpectraDataMod(FolderLoop,FileLoop).LockInNNTemp{1,iNNs}));
273	        end
274	
275	        %PHASE
276	        SpectraDataMod(FolderLoop,FileLoop).PhaseRef = ...
277	            angle(SpectraDataMod(1,1).fftRefGrn(ModFreqIndex));
278	        SpectraDataMod(FolderLoop,FileLoop).PhasePeak = ...
279	            angle(SpectraDataMod(FolderLoop,FileLoop).LockInPeak);
280	        SpectraDataMod(FolderLoop,FileLoop).PhaseIntegrated = ...
281	            angle(SpectraDataMod(FolderLoop,FileLoop).LockInIntegrated);
282	        SpectraDataMod(FolderLoop,FileLoop).PhasePeakWL = ...
283	            angle(SpectraDataMod(FolderLoop,FileLoop).LockInPeakWL);
284	        for iNNs = 1:NNs
285	            SpectraDataMod(FolderLoop,FileLoop).PhaseNNTemp{1,iNNs}  = ...
286	            angle(SpectraDataMod(FolderLoop,FileLoop).LockInNNTemp{1,iNNs});
287	        end
288	    end
289	end
290	
291	%% Put into file for fitmain that can be read
292	for FolderLoop = 1:length(Folder.TempCell)
293	    %Counts the number of 
294	    files = dir(Folder.TempCell{1,FolderLoop}); % equal
295	    files(strncmp({files.name}, '.', 1)) = []; % new, no exceptions
296	    
297	    for FileLoop = 1:length(files)
298	        %Position in m
299	        z(FileLoop,1)           = SpectraDataMod(FolderLoop,FileLoop).MotorPosition/1000;
300	        
301	        omega(FileLoop,1)       = SpectraDataMod(FolderLoop,FileLoop).ModFreq;
302	        
303	        MagPeak(FileLoop,1)     = SpectraDataMod(FolderLoop,FileLoop).LnMagPeak;
304	        MagInteg(FileLoop,1)    = SpectraDataMod(FolderLoop,FileLoop).LnMagIntegrated;
305	        MagPeakWL(FileLoop,1)   = SpectraDataMod(FolderLoop,FileLoop).LnMagPeakWL;
306	        for iNNs = 1:NNs
307	            MagNNTemp(FileLoop,iNNs) = SpectraDataMod(FolderLoop,FileLoop).LnMagNNTemp{1,iNNs};
308	        end
309	        
310	        PhaseRef(FileLoop,1)    = SpectraDataMod(FolderLoop,FileLoop).PhaseRef;
311	        ref = PhaseRef(FileLoop,1);
312	        PhasePeak(FileLoop,1)   = ref - SpectraDataMod(FolderLoop,FileLoop).PhasePeak;
313	        PhaseInteg(FileLoop,1)  = ref - SpectraDataMod(FolderLoop,FileLoop).PhaseIntegrated;
314	        PhasePeakWL(FileLoop,1) = ref - SpectraDataMod(FolderLoop,FileLoop).PhasePeakWL;
315	        for iNNs = 1:NNs
316	            PhaseNNTemp(FileLoop,iNNs) = ref - ...
317	            SpectraDataMod(FolderLoop,FileLoop).PhaseNNTemp{1,iNNs};
318	        end
319	    end
320	    
321	    cd(Folder.ForFitMain)
322	    TempDate=[SpectraDataMod(FolderLoop,FileLoop).Temperature,'K - ',FileData.monthStr,FileData.dateStr];
323	    MonthDate = [FileData.monthStr,FileData.dateStr];
324	    
325	    mkdir(MonthDate)
326	    cd([Folder.ForFitMain,'\',MonthDate])
327	    mkdir(TempDate)
328	    cd([Folder.ForFitMain,'\',MonthDate,'\',TempDate])
329	    
330	    [b,i] = sort(z);
331	    oldz = z;
332	    z=z(i);
333	        
334	    % write position and Magnitude into files for each temperature
335	    MagPeak=MagPeak(i);
336	    MagIngeg=MagInteg(i);
337	    MagPeakWL=MagPeakWL(i);
338	    dlmwrite([TempDate,', LnMagPeak.txt'], [z,MagPeak])
339	    dlmwrite([TempDate,', LnMagInteg.txt'], [z,MagInteg])
340	    dlmwrite([TempDate,', LnMagPeakWL.txt'], [z,MagPeakWL])
341	    for iNNs = 1:NNs
342	        MagNNTemp(:,iNNs)=MagNNTemp(i,iNNs);
343	        dlmwrite([TempDate,', LnMagNNTemp',num2str(iNNs),'.txt'], [z,MagNNTemp(i,iNNs)])
344	    end
345	    
346	    % write position and Phase into files
347	    PhasePeak=PhasePeak(i);
348	    PhaseIngeg=PhaseInteg(i);
349	    PhasePeakWL=PhasePeakWL(i);
350	    dlmwrite([TempDate,', PhasePeak.txt'], [z,PhasePeak])
351	    dlmwrite([TempDate,', PhaseInteg.txt'], [z,PhaseInteg])
352	    dlmwrite([TempDate,', PhasePeakWL.txt'], [z,PhasePeakWL])
353	    for iNNs = 1:NNs
354	        PhaseNNTemp(:,iNNs)=PhaseNNTemp(i,iNNs);
355	        dlmwrite([TempDate,', PhaseNNTemp',num2str(iNNs),'.txt'], [z,PhaseNNTemp(i,iNNs)])
356	    end
357	    
358	    % write position and frequencies
359	    dlmwrite([TempDate,', Frequencies.txt'], [z,omega])
360	end
361	
362	cd(Folder.SpectraExpHome)
363	filename = [FileName.DateStr,'.mat'];
364	save(filename,'-v7.3')
365	aveFreq = mean(omega);
366	close all
367	
368	
369	
370	
371	
372	%% Call Fitmain to calculate slope of MAG then PHASE
373	% Make thing to be saved that goes, mag of peak, integ, WL then phase, then
374	% temperature, then motor position, then frequency
375	cd(Folder.ForFitMain)
376	filename
377	beep;beep;beep;beep;beep;
378	[y,Fs] = audioread([Folder.Sound,'\Ready for Fitting.mp3']);sound(y,Fs);
379	pause
380	clearvars -except omega filename
381	
382	%%
383	params = [1600      1.75    0.0004      11];
384	fitmain
385	
386	%%
387	mLnMag = params(1);
388	aveFreq = mean(omega);
389	aveFreq=f;
390	
391	alpha = abs((pi*aveFreq)/mLnMag^2)
392	