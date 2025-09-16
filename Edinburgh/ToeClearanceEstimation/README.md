# Log of Toe Clearance Estimation Model


## 10/Jul/2025 - 14/Jul/2025
1. Some information missing in the `.mat` dataset, such as `hs_r_time`，`shokac_starting_time`, `time_diff`, and `Timesteps`, which are necessary for model training. So, a script, `DataExtraProcessing_xwt.m`, was write to supplement these missing information into the old dataset. But, note that there are two files in folder `S1-->ShoeData` being misspelled 'in' as 'n'. As a result, the script prompts errors when processing the two files. Just revise the file name and the script will go smoothly.

2. Some changes were made to the script `ReadingShoeDataTrainingModels_Xiaowei.m`, which help the convergence of model, but fine-tuning is still required. The next step is to try some different model parameters. There are some notes about the changes made to the script (not the all, see the new script for all changes):

  - 2.1 
    ```matlab
    shokac_time = results.(subjects{i}).ShokacData.(experiment_name).time(k+window) - ...
                  results.(subjects{i}).ShokacData.(experiment_name).time(1);
    ```
    the `time(k+window)` should be `time(k : k+window)`.

  - 2.2
    ```matlab
    eq_vicon_time = shokac_time + results.(subjects{i}).ShokacData.(experiment_name).time_diff;
    ```
    the `+` sign may should be `-` sign.

  - 2.3  
    If there is a frame of data that contain `NaN` elements, the whole frame should be deleted.

  - 2.4  
    Some parameters and structure of the model were modified to improve convergence and accuracy, but fine-tuning is still required.

  - 2.5  
    I don't know the meaning of `masked_data` and cannot parse it. Only S1 dataset has the `masked_data`. This information cannot be added into the other subjects' dataset via the script `DataExtraProcessing_xwt.m`. Currently, I just removed it from the model training and testing datset.

3. These may be some information that is more directly related to the toe clearance, such as the placement position of IMU, foot length, acc, gyro, and foot orientation angles. These information may help improve the model accuracy, but some of them were not provided by the instrumented shoes. I will try to obtain these information from Marker trajectory, as a test to know which types of info are most relevant or important for toe clearance estimation.


## 14/Jul/2025 - 22/Jul/2025
1. I guessed the most relevant info to estimate toe clearance includes: (1) 3D foot orientation angle, (2) 3D linear acceleration, (3) absoluate time (for integral), (4) foot length and the position of IMU on the foot, and (5) treadmill inclination. Thinking that the treadmill speed or walking on a treadmill only has an impact on the foot position estimation on the treadmill surface plane (left-right, forward-backward), which means that the toe clearance will not be affected by the treadmill speed itself.

2. I checked the time alignment between GRF and shoe data. They are mostly aligned, but the shoe data are a little instable.

3. I added the stance & swing phase informaton into the input dataset, which however did not improve the training results. Using only phase information gave worse results compared to using only pressure data, which might be right since stance phase does not mean that there will not be a change in the toe clearance value.

4. Toe clearance may should be the Marker `MTP1 Y-axis` value, but the toe clearance in the old dataset do not match the `MTP1_Y` value. Have no idea which one is incorrect. Also, the old dataset sometimes contains negative toe clearance values, which does not make sense.  
I tried using the `MTP1_Y` value (after subtracting the minimum to avoid negative elements) as the true toe clearance for training, but the training results become worser.  
   **Solved**: The reference toe clearance is not exactly the `MTP1 Y-axis` value, since the shoe shape has to be considered to obtain exactly the toe cleance. The marker `MTP1` measures the shoe clearance, not the toe segment. So, it was a hard work to calculate the reference toe clearance.

5. How is the toe clearance defined when the treadmill has a inclination? Is the perpendicular distance between the toe point against the treadmill surface? In the old dataset, toe clearance values increase during the stance phase, which seems incorrect. It should likely decrease.    
   **Solved**: Yes. It is the perpendicular distance.

6. The model training has a relative better results when only using dataset collected under zero inclination conditions (might indicate the issues in the toe clearance calculation for inclined surfaces).  
   **Solved**: No. Not the reference toe clearance calculation's fault. It is possibly because the lack of gyro or angle information in the dimensions expect sagittal plane, so as the measured body-frame-based acc data are not suitable anymore to estimate the abosulte trajectories.

7. I have tried different scale factors and model parameters. Up to now, the best results (both training and final test errors are around 10 mm) appear when using the combination below (using only zero-inclination dataset): (1) unnormalized time, (2) unnormalized acc in unit of m/s2, (3) unnormalized gyro in unit of deg/s, (4) unnormalized treadmill inclination in unit of deg, (5) unnormalized shoe size in unit of m, and (6) normalized foot pressure.

8. Questions:
  - 8.1 Can 3D gyro or, more directly, 3D orientation angle be achieved via the Delsys sensors?  
    **Solved**: Yes, but we want to use the data from the shoe only, so as to buid a portable system in the future. But, we can get other info, such as angles, from Marker data to test the potential optimal combination of model input data.
  
  - 8.2 How are the reference toe clearance values calculated? especially when there is a treadmill inclination.  
    **Solved**: As replied before. The calculation is correct.

9. Should the grf data from shoes have a very high accuracy? I mean if it is necessary to reduce the offset of the grf data. Actually, in physics, I don't think the pressure data play a important role in toe clearance estimation.  
   **Solved**: Maybe more accurate and more detailed grf data can help the model become better, but we dont know. Actually, the grf data are mainly used to separate the stance and swing phases.

10. using logic value for pressure data, to verify that the pressure data value dont need to be so accurate. I know that there are some offset in the collected grf data but that may have no effects on the estimation accuracy. (Threshld is set at 0.5 to convert force to logic values; this pipeline shows better results)

11. filter tc (3rd order butterworth, Fs 100Hz, Fc 15Hz; filtered data lead to better results)

12. Verify the correctness of the baseline toe clearance.  
    **Solved**: The current baseline toe clearances are correct, the difference in lowest value from different subjects is due to the deformation of the shoe becuase of different weight of subjects.

13. 3D angle data are necessary, especially when the acc data are presented in the body coordiante system; otherwise, the 3D angle data are not required here. (at least that there are 3D angle or grf data to help distinguish the stance or swing phase so as to know when the toe clearance is zero.)  
    **Solved**: Not sure the coordinate system that the acc data are represented in. We don't have 3D angle data currently from shoe. just a try of creating a new model, whose inputs are shoe data and outputs are the 3D angle data culculated from marker data. double models.

14. The gyro data from different axes have the same value. It is the fault of the sensor itself. But the 3D gyro data, or angle data, are pretty necessary to transform the acc data from the body-frame to the global-frame.



**Next Steps:**  
1. try to reduce offset in grf data to make them more accurate, so we expect that the model can learn some info from the accurate grf.
2. try to confirm the coordinate system that the acc data are represented in, body-based or world-based.
3. try to parse foot 3D angles from marker data and regard them as input data to test the accuracy. And then, trying to build a new model to predict 3D angles (not sure will have a good estimation). Maybe when we have 3D angles, the input dataset of collected in inclination will also produce a good results of training.
4. read the *Scientific Report* paper.


## 22/Jul/2025 - 29/Jul/2025
1. There is a trend that z axis acc value decreased and the y axis acc value increased when the inclination changed from zero to a non-zero value. And from the perspective of IMU itself, it essentially measures acc in the body-frame. If acc in world-frame is request, there must be an in-built algorithm embeded in the IMU sensor to transform the coordinate systems for acc data and this transformation must rely on 3D angles (if it is, the IMU should be able to output 3D angles).  
   **So, I think the acc measurement should be represented in the body-based coordinate system. Thereby, 3D angles are still required as input to the model so as to better estimate toe clearance.**

2. **The IMU 3D orientation angles might have been accurately estimated based on Marker data, as provided in** `calc3DIMUOrientation.m`.

3. Adding the 3D orientation angles considerably improve the model outcome. Moreover, including the data collected with non-zero treadmill inclication also help improves the outcome possibly due to the increased amount of training data. The results are given below:  
 
	   Train RMSE = 3.2271 mm  
	   Train R2 = 0.9653  
	   Train Bias = -0.05 mm  
	   Train 95CI = [-6.37 mm, 6.28 mm]  
	   Validation RMSE = 5.8550 mm  
	   Validation R2 = 0.8731  
	   Validation Bias = 0.73 mm  
	   Validation 95CI = [-10.65 mm, 12.12 mm]


4. Both side of data were used, so the amount of dataset has doubled, results (the same configuration with above) are as follows:
 
	    Train RMSE = 4.0920 mm  
	    Train R2 = 0.9460  
	    Train Bias = -0.20 mm  
	    Train LoA = [-8.21 mm, 7.81 mm]  
	    Validation RMSE = 6.2310 mm  
	    Validation R2 = 0.8737  
	    Validation Bias = -0.25 mm  
	    Validation LoA = [-12.45 mm, 11.96 mm]  



## 29/Jul/2025 - 29/Aug/2025
1. Adding the edge case data makes the model outcome worse (S1 validation; S2-S10 train and test), as follows:  

	 	Train RMSE = 4.4024 mm  
	    Train R2 = 0.9442  
	    Train Bias = -0.04 mm  
	    Train LoA = [-8.67 mm, 8.59 mm]  
	    Validation RMSE = 7.5906 mm  
	    Validation R2 = 0.7867  
	    Validation Bias = 1.08 mm  
	    Validation LoA = [-13.64 mm, 15.81 mm]


2. Using the calibrated shoe grf (zero the forces during swing phase) as the model input (real value, not logic value; without edge case data), obtains the results as follows:
 
	    Train RMSE = 3.0697 mm  
	    Train R2 = 0.9686  
	    Train Bias = -0.07 mm  
	    Train 95CI = [-6.09 mm, 5.94 mm]  
	    Validation RMSE = 5.5402 mm  
	    Validation R2 = 0.8862  
	    Validation Bias = 1.02 mm  
	    Validation 95CI = [-9.65 mm, 11.70 mm]  


3. the same configuration with the above, only the calibrated shoe grfs are converted into logic values (0.5N or 0.5Nm threshold):  

	 	Train RMSE = 3.1349 mm  
	    Train R2 = 0.9672  
	    Train Bias = -0.02 mm  
	    Train 95CI = [-6.16 mm, 6.12 mm]  
	    Validation RMSE = 5.6893 mm  
	    Validation R2 = 0.8800  
	    Validation Bias = 1.12 mm  
	    Validation 95CI = [-9.82 mm, 12.05 mm]  


4. Leave-one-subject-out results of each subject as the validation and the remaining nine subjects as the training and test (without including edge case data; using the new version of code [01-Aug-2025]):

|Validation Subject|S1|S2|S3|S4|S5|S6|S7|S8|S9|S10|
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|Train RMSE /mm|2.3554|2.2950|2.5872|2.6026|2.7469|2.3812|2.1220|2.2416|2.7126|2.2487|
|Train R2|0.9815|0.9821|0.9777|0.9775|0.9755|0.9810|0.9852|0.9826|0.9756|0.9822|
|Train Bias /mm|-0.06|0.00|-0.13|0.00|-0.17|0.01|0.04|0.07|-0.03|0.00|
|Train 95CI /mm|[-4.67, 4.56]|[-4.50, 4.49]|[-5.19, 4.94]|[-5.10, 5.11]|[-5.55, 5.20]|[-4.66, 4.67]|[-4.12, 4.20]|[-4.33, 4.46]|[-5.34, 5.29]|[-4.41, 4.41]|
|Validation RMSE /mm|4.8271|3.9953|6.0547|3.7665|4.9172|5.7753|3.6737|4.9109|4.4955|5.5455|
|Validation R2|0.9136|0.9495|0.8667|0.9483|0.8754|0.8858|0.9461|0.9381|0.9246|0.9240|
|Validation Bias /mm|0.52|-1.94|2.44|-0.88|-1.03|1.29|-1.06|-0.90|2.46|1.98|
|Validation 95CI /mm|[-8.88, 9.93]|[-8.79, 4.90]|[-8.43, 13.30]|[-8.06, 6.29]|[-10.46, 8.39]|[-9.74, 12.32]|[-7.95, 5.84]|[-10.36, 8.56]|[-4.91, 9.83]|[-8.17, 12.13]|
|Training Duration /min|263|277|174|222|220|264|410|286|203|333|

The foot IMU angles used in the above LOSO test are incorrect due to the uses of absoluate index values in marker column data extraction and also frames that were inproperly defined. Frames involved are redefined as presented in file 'Frames' in this repo and angles are re-estimated (have been double checked). Below are the results of LOSO retest without any changes to the script except using the correct foot angles (note the both LOSO tests used the non-boolean GRFs data, 0.02 scale_y, and bilstm, mentioned in the next item). All results are improved, and some changes are significant.

|Validation Subject|S1|S2|S3|S4|S5|S6|S7|S8|S9|S10|
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|Train RMSE /mm			|1.8878|2.2647|2.0154|2.2316|2.5319|1.9511|2.1039|2.0658|2.5412|2.4696|
|Train R2				|0.9881|0.9826|0.9864|0.9834|0.9792|0.9873|0.9854|0.9852|0.9786|0.9785|
|Train Bias /mm			|-0.02|-0.16|-0.07|-0.02|0.06|0.01|-0.00|0.04|-0.01|-0.13|
|Train 95CI /mm			|[-3.72, 3.68]|[-4.59, 4.27]|[-4.02, 3.88]|[-4.39, 4.36]|[-4.90, 5.03]|[-3.82, 3.83]|[-4.13, 4.12]|[-4.01, 4.09]|[-4.99, 4.97]|[-4.97, 4.70]|
|Validation RMSE /mm	|4.2389|3.5448|4.5685|3.1391|3.8198|3.9789|3.3390|4.7764|3.4994|4.5941|
|Validation R2			|0.9334|0.9603|0.9241|0.9641|0.9248|0.9458|0.9554|0.9414|0.9543|0.9478|
|Validation Bias /mm	|-0.27|-1.56|2.07|-0.97|-0.48|1.13|-1.09|-0.96|1.80|1.94|
|Validation 95CI /mm	|[-8.56, 8.02]|[-7.80, 4.68]|[-5.91, 10.05]|[-6.82, 4.89]|[-7.91, 6.94]|[-6.35, 8.61]|[-7.28, 5.10]|[-10.13, 8.21]|[-4.08, 7.68]|[-6.23, 10.10]|
|Training Duration /min	| | | | | | | | | | |

Please note adding `'ExecutionEnvironment', 'auto'` in the `trainingOptions` can reduce the training duration by more than half.  

If adding the shank angle into the input. Seem get the results worse.   

|Validation Subject|S1|S2|S3|S4|S5|S6|S7|S8|S9|S10|
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|Train RMSE /mm			|2.1685| | | | | | | | |2.1951|
|Train R2				|0.9843| | | | | | | | |0.9830|
|Train Bias /mm			|-0.10| | | | | | | | |0.02|
|Train 95CI /mm			|[-4.35, 4.15]| | | | | | | | |[-4.29, 4.32]|
|Validation RMSE /mm	|4.9081| | | | | | | | |5.5249|
|Validation R2			|0.9107| | | | | | | | |0.9245|
|Validation Bias /mm	|-1.80| | | | | | | | |3.16|
|Validation 95CI /mm	|[-10.75, 7.16]| | | | | | | | |[-5.72, 12.04]|
|Training Duration /min	| | | | | | | | | | |


5. Changing the scale_y to 0.02, adding the magnitude of the Shokac GRFs in the training (not just the booleans), and making the first layer of the network a 'bilstm' layer:

		Train RMSE = 2.2845 mm  
		Train R2 = ...  
		Train Bias = ...  
		Train 95CI = ...  
		Validation RMSE = 4.6970 mm  
		Validation R2 = 0.9182  
		Validation Bias = 0.58 mm  
		Validation 95CI = [-8.55 mm, 9.72 mm]  

6. Test edge case data:  
   Validation: S1; Train: the first 80% of S2-S10; Test: the last 20%.  
   Edge Train: 20%-80% of all edge case data; Edge Test: the last 20%; Edge Validation: the first 20%.  

   Train the model using normal data first, and then fine tune the model obtained in the first step using the Edge Train and Edge Test data.  
   
Pipeline #1: train the second model without freezing any parameters or layers of the first model, but just decrease the learning rate to 1/50 of its value used in the first training process. The before and after below mean that before and after the fine tuning:

|Types|Train RMSE|Train R2|Train Bias|Train 95CI|Val RMSE|Val R2|Val Bias|Val 95CI|Edge RMSE|Edge R2|Edge Bias|Edge 95CI|
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
|Before|2.3554|0.9815|-0.06|[-4.67, 4.56]|4.8271|0.9136|0.52|[-8.88, 9.93]|34.3661|0.5552|-9.12|[-74.06, 55.83]|
|After|19.1601|-0.2235|6.85|[-28.22, 41.92]|12.5774|0.4137|3.18|[-20.67, 27.03]|34.9673|0.5395|-8.64|[-75.05, 57.77]|

Pipeline #2: train the second model while freezing the first three layers of the first model, any parameters and weights in the first three layers are allowed to be changed (have tested without any freezes to the first model, but achieved considerably worse results, since the first model has been disrupted by the edge case data). the learning rate used in the second training is 1/10 of the previous one.

|Types|Train RMSE|Train R2|Train Bias|Train 95CI|Val RMSE|Val R2|Val Bias|Val 95CI|Edge RMSE|Edge R2|Edge Bias|Edge 95CI|
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
|Before|2.3554|0.9815|-0.06|[-4.67, 4.56]|4.8271|0.9136|0.52|[-8.88, 9.93]|34.3661|0.5552|-9.12|[-74.06, 55.83]|
|After|50.7679|-7.5896|15.62|[-79.05, 110.30]|25.8999|-1.4862|6.03|[-43.34, 55.40]|37.3556|0.4744|-16.94|[-82.20, 48.32]|


7. Test results under various configurations

|Condition|A|A|A|A|B|C|C|D|
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
|V-Subject	|S1|S3|S5|S7|S1|S1|S10|S1|
|T-RMSE /mm	|2.2411|2.0813|2.0372|2.5285|2.1984|2.3289|2.6756|2.5171|
|T-R2		|0.9833|0.9855|0.9865|0.9789|0.9839|0.9819|0.9748|0.9789|
|T-Bias /mm	|0.02|-0.08|-0.13|-0.14|-0.00|-0.07|-0.03|0.16|
|T-95CI /mm	|[-4.37, 4.42]|[-4.15, 4.00]|[-4.11, 3.86]|[-5.09, 4.81]|[-4.31, 4.31]|[-4.63, 4.49]|[-5.27, 5.22]|[-4.76, 5.08]|
|V-RMSE /mm	|4.5559|4.7148|3.6599|3.7455|4.4249|5.0858|4.9275|4.1496|
|V-R2		|0.9231|0.9192|0.9310|0.9439|0.9274|0.9041|0.9400|0.9362|
|V-Bias /mm	|-0.63|2.33|-0.56|-1.15|-0.10|-0.41|1.62|0.35|
|V-95CI /mm	|[-9.47, 8.21]|[-5.71, 10.36]|[-7.65, 6.53]|[-8.14, 5.83]|[-8.77, 8.58]|[-10.35, 9.53]|[-7.50, 10.74]|[-7.75, 8.46]|

	- A: added masked_data
	- B: added logical phase index (1-stance, 0-swing)
 	- C: normalized time, shoe grf, acc, gyro, inc, shoe size
  	- D: normalized shoe grf, shoe size


|Condition|A|B|C|D|E|F|G|
|:---|---:|---:|---:|---:|---:|---:|---:|
|V-Subject	|S1|S1|S1|S1|S1|S1|S1|
|T-RMSE /mm	|8.9158|2.2182|2.1796|1.9381|2.3729|1.9405|1.7826|
|T-R2		|0.7351|0.9836|0.9842|0.9875|0.9812|0.9875|0.9894|
|T-Bias /mm	|0.24|-0.18|-0.01|-0.02|0.11|-0.03|0.01|
|T-95CI /mm	|[-17.23, 17.71]|[-4.51, 4.15]|[-4.28, 4.27]|[-3.82, 3.78]| [-4.53, 4.76]|[-3.83, 3.77]|[-3.48, 3.51]|
|V-RMSE /mm	|7.9720|4.3832|4.4779|4.2451|4.4801|4.3871|4.7259|
|V-R2		|0.7645|0.9288|0.9257|0.9332|0.9256|0.9287|0.9172|
|V-Bias /mm	|1.02|-0.63|-0.29|-0.22|-0.06|-0.40|-0.96|
|V-95CI /mm	|[-14.48, 16.52]|[-9.13, 7.87]|[-9.05, 8.47]|[-8.53, 8.09]|[-8.84, 8.72]|[-8.96, 8.16]|[-10.03, 8.11]|

	- A: added masked_data, added body weight, changed shoe length to 250mm (40) and 263mm (42), initial value of time series wasn't subtracted
	- B: A - time series changes
 	- C: B - body weight
  	- D: C - masked_data + normalize time series (0-1)
    - E: C - masked_data
	- F: changed shoe size from 40 and 42 to 285mm and 300mm (measured from STL files)
 	- G: changed `step_between_points` from `window/10` to `round(window/8)`

 |Condition|A|B|C|D|E|F| |
|:---|---:|---:|---:|---:|---:|---:|---:|
|V-Subject	|S1|S1|S1|S1|S1|S1| |
|T-RMSE /mm	|2.5629|1.8878|2.2829|2.4184|1.8542|2.3069| |
|T-R2		|0.9781|0.9881|0.9826|0.9805|0.9885|0.9823| |
|T-Bias /mm	|-0.08|-0.02|-0.07|0.28|-0.08|0.03| |
|T-95CI /mm	|[-5.10, 4.94]|[-3.72, 3.68]|[-4.54, 4.40]|[-4.43, 4.99]|[-3.71, 3.55]|[-4.49, 4.55]| |
|V-RMSE /mm	|4.8336|4.2389|4.5388|4.9301|4.3035|4.6680| |
|V-R2		|0.9134|0.9334|0.9236|0.9099|0.9314|0.9192| |
|V-Bias /mm	|-0.55|-0.27|-0.17|-0.34|-0.32|-0.10| |
|V-95CI /mm	|[-9.97, 8.86]|[-8.56, 8.02]|[-9.06, 8.72]|[-9.98, 9.30]|[-8.73, 8.09]|[-9.24, 9.05]| |

	- A: changed `L2Regularization` from 1e-4 to 9e-4 (to mitigate over-fitting)
	- B: changed `L2Regularization` to 5e-4
 	- C: B + changed `dropoutLayer1` from 0.2 to 0.3
  	- D: `dropoutLayer1` 0.25, `L2Regularization` 2e-4 
   	- E: numHiddenUnits1 = numFeatures*4; numHiddenUnits2 = numFeatures*2; mini_batch_size=256; 'ValidationPatience', 15;
	- F: removed logical shoe GRF values

 |Condition|A|A|A|A|A|B|C|D|
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
|V-Subject	|S1|S3|S5|S7|S9|S1|S1|0%-10%|
|T-RMSE /mm	|1.7358|1.9225|2.0063|2.0070|1.5938|1.9459|4.3397|1.8910|
|T-R2		|0.9900|0.9877|0.9869|0.9867|0.9916|0.9874|0.9372|0.9875|
|T-Bias /mm	|0.09|-0.08|0.14|0.15|0.20|-0.28|-0.44|-0.13|
|T-95CI /mm	|[-3.31, 3.49]|[-3.85, 3.68]|[-3.78, 4.06]|[-3.77, 4.07]|[-2.9, 3.3]|[-4.06, 3.49]|[-8.9, 8.02]|[-3.83, 3.5]|
|V-RMSE /mm	|3.7652|4.9587|3.0674|3.0935|3.5234|3.9014|6.1124|4.2867|
|V-R2		|0.9475|0.9106|0.9515|0.9618|0.9537|0.9436|0.8615|0.9308|
|V-Bias /mm	|0.05|2.02|-0.37|-0.52|2.04|-0.51|-0.37|0.11|
|V-95CI /mm	|[-7.33, 7.43]|[-6.86, 10.9]|[-6.34, 5.6]|[-6.5, 5.45]|[-3.6, 7.67]|[-8.09, 7.07]|[-12.33, 11.59]|[-8.29, 8.51]|

	- A: changed model parameters as follows (added GELU):
  ```matlab
% changed to two-layer network, added geluLayer.
layers = [
    sequenceInputLayer(numFeatures) % input layer
    bilstmLayer(numFeatures*6, 'OutputMode', 'sequence')  % 1st layer
    batchNormalizationLayer
    geluLayer
    dropoutLayer(0.15)

    lstmLayer(numFeatures*3, 'OutputMode', 'sequence')  % 2nd layer
    batchNormalizationLayer
    geluLayer
    dropoutLayer(0.15)

    fullyConnectedLayer(numResponses) % output layer
    regressionLayer
    ];

options = trainingOptions('adam',...
    'ValidationData',{x_test, y_test},...
    'ValidationPatience', 12,... % changed
    'MiniBatchSize', mini_batch_size,...
    'MaxEpochs', NEpoch, ...
    'InitialLearnRate', LearningRate, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 4, ... % changed
    'LearnRateDropFactor', 0.6, ... % changed
    'GradientThreshold', 1, ... % added
    'L2Regularization', 1e-4, ...
    'Shuffle', shuffle, ...
    'Verbose', false, ...
    'Plots', 'training-progress',...
    'ExecutionEnvironment', 'auto');
  ```

	- B: changed model parameters as follows (out1 = bilstm; out2 = lstm; out2 = out1 + out2):
 ```matlab
layers = [
    sequenceInputLayer(numFeatures, 'Name', 'input') % Input layer
    bilstmLayer(numFeatures*6, 'OutputMode', 'sequence', 'Name', 'bilstm') % 1st layer
    batchNormalizationLayer('Name', 'bn1')
    geluLayer('Name', 'gelu1')
    dropoutLayer(0.15, 'Name', 'dropout1')

    lstmLayer(numFeatures*3, 'OutputMode', 'sequence', 'Name', 'lstm') % 2nd layer
    batchNormalizationLayer('Name', 'bn2')
    geluLayer('Name', 'gelu2')
    dropoutLayer(0.15, 'Name', 'dropout2')
    fullyConnectedLayer(numFeatures*12, 'Name', 'fc_lstm') % Adjust dimension to match BiLSTM

    additionLayer(2, 'Name', 'add') % Residual connection (two inputs), place here to add the output of above, i.e., out2
    batchNormalizationLayer('Name', 'bn3')
    geluLayer('Name', 'gelu3')

    fullyConnectedLayer(numResponses, 'Name', 'fc_output') % Output layer
    regressionLayer('Name', 'output')
];

% Create layer graph
lgraph = layerGraph(layers);
% Add residual connection from bilstm to additionLayer
lgraph = connectLayers(lgraph, 'bilstm', 'add/in2'); % to add out1

% Training options (unchanged)
options = trainingOptions('adam', ...
    'ValidationData', {x_test, y_test}, ...
    'ValidationPatience', 12, ...
    'MiniBatchSize', mini_batch_size, ...
    'MaxEpochs', NEpoch, ...
    'InitialLearnRate', LearningRate, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 4, ...
    'LearnRateDropFactor', 0.6, ...
    'GradientThreshold', 1, ...
    'L2Regularization', 1e-4, ...
    'Shuffle', shuffle, ...
    'Verbose', false, ...
    'Plots', 'training-progress', ...
    'ExecutionEnvironment', 'auto');

% Train the network
[net, info] = trainNetwork(x_train, y_train, lgraph, options);
```

 	- C: changed model parameters as follows (out1 = bilstm; out2 = out1 + inputs):
 ```matlab
layers = [
    sequenceInputLayer(numFeatures, 'Name', 'input') % Input layer
    fullyConnectedLayer(numFeatures*12, 'Name', 'fc_input') % Adjust input dimension to match BiLSTM
    
    bilstmLayer(numFeatures*6, 'OutputMode', 'sequence', 'Name', 'bilstm') % BiLSTM layer
    
    additionLayer(2, 'Name', 'add') % Residual connection (two inputs), place here to add the output of above, i.e., out1
    batchNormalizationLayer('Name', 'bn1')
    geluLayer('Name', 'gelu1')
    dropoutLayer(0.15, 'Name', 'dropout1')

    lstmLayer(numFeatures*3, 'OutputMode', 'sequence', 'Name', 'lstm') % LSTM layer
    batchNormalizationLayer('Name', 'bn2')
    geluLayer('Name', 'gelu2')
    dropoutLayer(0.15, 'Name', 'dropout2')

    fullyConnectedLayer(numResponses, 'Name', 'fc_output') % Output layer
    regressionLayer('Name', 'output')
];

% Create layer graph
lgraph = layerGraph(layers);
% Add residual connection from input to additionLayer
lgraph = connectLayers(lgraph, 'fc_input', 'add/in2'); % add the modified input that has the same size structure

% Training options (unchanged)
options = trainingOptions('adam', ...
    'ValidationData', {x_test, y_test}, ...
    'ValidationPatience', 12, ...
    'MiniBatchSize', mini_batch_size, ...
    'MaxEpochs', NEpoch, ...
    'InitialLearnRate', LearningRate, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 4, ...
    'LearnRateDropFactor', 0.6, ...
    'GradientThreshold', 1, ...
    'L2Regularization', 1e-4, ...
    'Shuffle', shuffle, ...
    'Verbose', false, ...
    'Plots', 'training-progress', ...
    'ExecutionEnvironment', 'auto');

% Train the network
[net, info] = trainNetwork(x_train, y_train, lgraph, options);
```
 	- D: model parameters are same as in Condition A, while all input data are z-score normalized:
  		- (1) collect each type of input/output data across all subjects and trials into a vector to calculate overall `mean` and `std`
		- (2) z-score normalize each type of input/output data, e.g., `(tc_l - mean_tc_l) / std_tc_l`
  		- (3) divide normalized data into moving winodws, and allocate 0%-10% for validation, 10%-80% for training, 80%-100% for testing
		- (4) codes afer `%% ↓ remove NaN elements, which may affects the model training` are unchanged


 |Condition|A|B|C|D|E|F|G|H|
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
|V-Subject	|S1|S1|S1|S1|S1|S1|S1|S1|
|T-RMSE /mm	|1.5963|1.8689|2.1802|2.3159|1.8189|1.7881|2.1339|2.5121|
|T-R2		|0.9915|0.9884|0.9842|0.9821|0.9890|0.9893|0.9848|0.9790|
|T-Bias /mm	|0.19|0.07|-0.12|-0.05|0.03|-0.43|0.03|-1.06|
|T-95CI /mm	|[-2.91, 3.30]|[-3.59, 3.73]|[-4.39, 4.15]|[-4.58, 4.49]|[-3.53, 3.59]|[-3.83, 2.97]|[-4.15, 4.21]|[-5.52, 3.41]|
|V-RMSE /mm	|3.9939|4.3660|4.2874|4.3910|3.9307|3.9748|4.1142|4.7835|
|V-R2		|0.9409|0.9293|0.9319|0.9285|0.9427|0.9414|0.9373|0.9152|
|V-Bias /mm	|0.02|-0.46|-0.05|-0.24|-0.25|-0.54|-0.45|-1.48|
|V-95CI /mm	|[-7.81, 7.85]|[-8.97, 8.05]|[-8.45, 8.36]|[-8.83, 8.36]|[-7.94, 7.44]|[-8.26, 7.18]|[-8.47, 7.56]|[-10.40, 7.44]|

 	- out1 = bilstm, out2=lstm; out2=out1+out2. It has the following combinations:
  	- A: bilstm → BN1 → GELU1 → Dropout1 → LSTM → +bilstm Output → regressionLayer
   	- B: bilstm → BN1 → GELU1 → Dropout1 → LSTM → +bilstm Output → BN2 → GELU2 → Dropout2 → regressionLayer
	- C: bilstm → BN1 → GELU1 → Dropout1 → LSTM → +Dropout1 Output → regressionLayer
 	- D: bilstm → BN1 → GELU1 → Dropout1 → LSTM → +Dropout1 Output → BN2 → GELU2 → Dropout2 → regressionLayer
  	- E: bilstm → BN1 → GELU1 → Dropout1 → LSTM → BN2 → GELU2 → Dropout2 → +bilstm Output → regressionLayer
    - F: bilstm → BN1 → GELU1 → Dropout1 → LSTM → BN2 → GELU2 → Dropout2 → +Dropout1 Output → regressionLayer
	- G: bilstm → LSTM → +bilstm Output → BN2 → GELU2 → Dropout2 → regressionLayer
 	- H: bilstm → LSTM → BN2 → GELU2 → Dropout2 → +bilstm Output → regressionLayer

  
-----
8. Test different input combinations

 |Condition|C1|C2|C3|C4|C5|C6|C7|C8|
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
|V-Subject	|0-10%|0-10%|0-10%|0-10%|0-10%|0-10%|0-10%|0-10%|
|T-RMSE /mm	|1.5082|1.6274|1.5889|1.9071|1.4121|1.4971|1.4055|1.9033|
|T-R2		|0.9920|0.9907|0.9912|0.9873|0.9930|0.9921|0.9931|0.9873|
|T-Bias /mm	|-0.05|-0.34|-0.20|0.02|-0.19|-0.05|0.13|-0.01|
|T-95CI /mm	|[-3.01, 2.90]|[-3.46, 2.78]|[-3.29, 2.89]|[-3.71, 3.76]|[-2.93, 2.55]|[-2.98, 2.88]|[-2.61, 2.87]|[-3.74, 3.72]|
|V-RMSE /mm	|4.3841|4.0704|4.2960|4.6544|4.1970|4.3340|4.2213|3.7260|
|V-R2		|0.9276|0.9376|0.9305|0.9184|0.9337|0.9293|0.9329|0.9477|
|V-Bias /mm	|-1.24|-1.38|-1.16|-1.05|-0.87|-1.06|-0.51|-0.24|
|V-95CI /mm	|[-9.48, 7.00]|[-8.89, 6.13]|[-9.27, 6.95]|[-9.94, 7.83]|[-8.92, 7.18]|[-9.29, 7.18]|[-8.72, 7.71]|[-7.52, 7.05]|


 |Condition|C9|C10|C11|C12|C13|C14|C15|C16|
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
|V-Subject	|0-10%|0-10%|0-10%|0-10%|0-10%|0-10%|0-10%|0-10%|
|T-RMSE /mm	|2.1535|1.4095|1.7660|2.0179|1.4199|1.8851|2.2805|1.4242|
|T-R2		|0.9838|0.9930|0.9891|0.9857|0.9929|0.9875|0.9818|0.9929|
|T-Bias /mm	|0.27|0.02|-0.46|0.53|-0.17|-0.29|0.25|-0.25|
|T-95CI /mm	|[-3.92, 4.45]|[-2.74, 2.78]|[-3.80, 2.88]|[-3.28, 4.35]|[-2.93, 2.59]|[-3.94, 3.36]|[-4.19, 4.70]|[-3.00, 2.50]|
|V-RMSE /mm	|4.1841|4.6216|4.6562|4.0231|4.2131|4.5442|5.1560|3.7065|
|V-R2		|0.9341|0.9196|0.9183|0.9390|0.9331|0.9222|0.8999|0.9483|
|V-Bias /mm	|0.55|-0.95|-1.15|0.78|-0.98|-1.09|1.01|-0.11|
|V-95CI /mm	|[-7.58, 8.68]|[-9.81, 7.92]|[-9.99, 7.70]|[-6.95, 8.52]|[-9.01, 7.05]|[-9.74, 7.56]|[-8.90, 10.92]|[-7.37, 7.15]|

 |Condition|C17|C18|C19|C20|C21|C22|C23|C24|
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
|V-Subject	|0-10%|0-10%|0-10%|0-10%|0-10%|0-10%|0-10%|0-10%|
|T-RMSE /mm	|1.6043|1.9783|1.7557|1.7429|2.1251|1.9671|1.5668|1.3376|
|T-R2		|0.9910|0.9863|0.9892|0.9894|0.9842|0.9864|0.9914|0.9937|
|T-Bias /mm	|-0.07|0.10|-0.27|-0.46|-0.36|-0.46|-0.17|-0.20|
|T-95CI /mm	|[-3.21, 3.07]|[-3.77, 3.97]|[-3.67, 3.13]|[-3.76, 2.83]|[-4.47, 3.74]|[-4.21, 3.29]|[-3.22, 2.88]|[-2.80, 2.39]|
|V-RMSE /mm	|4.4437|3.8905|4.5894|3.8079|4.4403|4.3107|4.1891|3.6511|
|V-R2		|0.9256|0.9430|0.9207|0.9454|0.9257|0.9300|0.9339|0.9498|
|V-Bias /mm	|-0.80|-0.59|-0.81|-0.91|-0.87|-1.23|-0.81|0.09|
|V-95CI /mm	|[-9.36, 7.77]|[-8.12, 6.95]|[-9.66, 8.05]|[-8.16, 6.34]|[-9.41, 7.66]|[-9.33, 6.87]|[-8.87, 7.25]|[-7.06, 7.25]|

 |Condition|C25|C26|C27| | | | | |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
|V-Subject	|0-10%|0-10%|0-10%| | | | | |
|T-RMSE /mm	|2.2663|1.8518|1.5804| | | | | |
|T-R2		|0.9820|0.9880|0.9912| | | | | |
|T-Bias /mm	|0.19|-0.16|-0.06| | | | | |
|T-95CI /mm	|[-4.23, 4.62]|[-3.78, 3.46]|[-3.15, 3.04]| | | | | |
|V-RMSE /mm	|3.7925|3.8035|3.3089| | | | | |
|V-R2		|0.9458|0.9455|0.9588| | | | | |
|V-Bias /mm	|0.48|0.45|0.46| | | | | |
|V-95CI /mm	|[-6.89, 7.85]|[-6.95, 7.86]|[-5.97, 6.88]| | | | | |


```matlab
% the completed inputs
nAll_Input.shoeTime(index, :)'; ...
nAll_Input.accel_l(index, :)'; ...
nAll_Input.EulerAngles_footIMU_l(index, :)';...
nAll_Input.EulerAngles_shank_l(index, :)';...
nAll_Input.gyro_l(index, :)';...
nAll_Input.heel_l(index, :)';...
nAll_Input.inc(index, :)';...
nAll_Input.little_l(index, :)';...
nAll_Input.shoeSize(index, :)';...
nAll_Input.thumb_l(index, :)']};

% model parameters usd
numFeatures = size(x_train{1}, 1);
numHiddenUnits1 = numFeatures*6;
numHiddenUnits2 = numFeatures*3;
numResponses = size(y_train{1}, 1);
dropoutLayer1=0.15;
NEpoch=80;
LearningRate=0.002;
shuffle = 'every-epoch'; % sequence-to-sequence training and estimation
layer1='sequence';
mini_batch_size=128;

layers = [
    sequenceInputLayer(numFeatures) % input layer
    bilstmLayer(numHiddenUnits1, 'OutputMode', 'sequence')  % 1st layer
    batchNormalizationLayer
    geluLayer
    dropoutLayer(dropoutLayer1)

    lstmLayer(numHiddenUnits2, 'OutputMode', 'sequence')  % 2nd layer
    batchNormalizationLayer
    geluLayer
    dropoutLayer(dropoutLayer1)

    fullyConnectedLayer(numResponses) % output layer
    regressionLayer
    ];

options = trainingOptions('adam',...
    'ValidationData',{x_test, y_test},...
    'ValidationPatience', 12,... % changed
    'MiniBatchSize', mini_batch_size,...
    'MaxEpochs', NEpoch, ...
    'InitialLearnRate', LearningRate, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 4, ... % changed
    'LearnRateDropFactor', 0.6, ... % changed
    'GradientThreshold', 1, ... % added
    'L2Regularization', 1e-4, ...
    'Shuffle', shuffle, ...
    'Verbose', false, ...
    'Plots', 'training-progress',...
    'ExecutionEnvironment', 'auto');

[net, info] = trainNetwork(x_train, y_train, layers, options);
```
 
	- C1: grf only has heel
 	- C2: grf only has little
  	- C3: grf only has thumb
    - C4: accel only has x,y
	- C5: accel only has x,z
 	- C6: accel only has y,z
  	- C7: foot angle only has x,y
    - C8: foot angle only has x,z
	- C9: foot angle only has y,z
 	- C10: shank angle only has x,y
 	- C11: shank angle only has x,z
 	- C12: shank angle only has y,z
  	- C13: remove time
    - C14: remove accel
	- C15: remove foot angle
 	- C16: remove shank angle
    - C17: remove gyro
	- C18: remove grf_heel
 	- C19: remove inc
  	- C20: remove grf_little
   	- C21: remove shoe size
	- C22: remove grf_thumb
 	- C23: completed inputs
  	- C24: foot angle only has x,z; remove shank angle
    - C25: foot angle only has x,z; remove shank angle; remove grf_heel
    - C26: foot angle only has x,z; remove shank angle; remove grf_little
    - C27: foot angle only has x,z; remove shank angle; remove grf_heel and grf_little
 
  
-----
9. Bayesian Optimisation for hyperparameter tuning:

```
	Best results:
		Train RMSE = 1.0186 mm
		Train R2 = 0.9964 
		Train Bias = -0.03 mm
		Train 95CI = [-2.02 mm, 1.97 mm]
		Validation RMSE = 3.8605 mm
		Validation R2 = 0.9439 
		Validation Bias = 0.03 mm
		Validation 95CI = [-7.53 mm, 7.60 mm]

	Search Space:
		numHiddenUnits1 = [round(0.25*numFeatures), round(0.5*numFeatures), numFeatures, 2*numFeatures, 4*numFeatures, 8*numFeatures];
		numHiddenUnits2 = [round(0.25*numFeatures), round(0.5*numFeatures), numFeatures, 2*numFeatures, 4*numFeatures, 8*numFeatures];
		dropout_rate = [0.10, 0.15, 0.2, 0.25, 0.30];
		LearningRate = [1e-05, 5e-05, 1e-4, 5e-04, 1e-03, 5e-03, 1e-2];
		L2Regularization =[1e-5, 1e-04, 1e-3];
		mini_batch_size = [64, 128, 256];

	Best model:
		numHiddenUnits1 = [8*numFeatures];
		numHiddenUnits2 = [8*numFeatures];
		dropout_rate = [0.10];
		LearningRate = [1e-2];
		L2Regularization =[1e-5];
		mini_batch_size = [256];
```
###
	Convergence figure (Normalised RMSE values):
<img width="3840" height="1868" alt="OptConvergenceFig" src="https://github.com/user-attachments/assets/182a0e41-fa58-4561-bab4-d99613f5d470" />


## 29/Aug/2025 - 08/Sep/2025
We conducted physical experiments to test the performance of our latest estimation model. However, the Normal Model (the model trained using normal data only) didn't perform well on extreme walking patterns. Edge case data must be utilized to address this issue.

1. Since edge case data aren't normal, they were z-score normalized using the mean and std paras calculated from the normal data.
2. Left and right data of edge cases weren't split to increase the size of available data for training. This change was also made in the Normal Model training process and a boolean variable indicating the body side was also added as one of the model inputs. This change doubled the data size and improved the Normal Model's performance.
3. Fine tunning was conducted using edge-normal mixed data, since training using only edge cases would cause catastrophic forgetting. The sequence of edge and normal data windows was shuffled before training to avoid biases. [tried to add the shuffle pipeline in the Normal Model training process and found it further improved the model performance, likely due to the enhanced generalizability to the sequence of subjects and conditions.]
4. Incorporating edge cases would inevitably degrade model performance on normal data. Sufficient edge and normal data may help find a optimal balance, but due to the limited dataset, a manual curriculum learning strategy was adopted to explore a model that has a acceptable performance on both edge and normal case data.  
   	4.1 Train the Normal Model seven times using datasets with 5%, 10%, 15%, 20%, 25%, 30%, and 35% edge case proportions. Each training was applied to the lastest model from the previous stage and used all the edge case data but the size of normal data was adjusted at each stage to achieve the desired edge-to-normal ratio.  
    4.2 The above pipeline was shown to considerably improve model performance compared to training once with a dataset containing 20% edge cases. Note that the normal data extracted to generate the mixed dataset at each stage were different. The sequence of edge and normal data windows was shuffled and different across stages.

```matlab
% fine tunning
NEpoch_ft = 20;
LearningRate_ft = 1e-3;
mini_batch_size_ft = 256;

options = trainingOptions('adam', ...
    'ValidationData',{x_val, y_val}, ...
    'ValidationPatience', 8, ...
    'MiniBatchSize', mini_batch_size_ft, ...
    'MaxEpochs', NEpoch_ft, ...
    'InitialLearnRate', LearningRate_ft, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropPeriod', 4, ...
    'LearnRateDropFactor', 0.6, ...
    'GradientThreshold', 1, ...
    'L2Regularization', 1e-3, ...
    'Shuffle','every-epoch', ...
    'Verbose',false, ...
    'Plots','training-progress', ...
    'ExecutionEnvironment','auto', ...
    'OutputNetwork','best-validation-loss');

if ~exist('net', 'var')
    load("net_base.mat"); % Normal Model
    layers_ft = net_base.Layers;
else
    net_last = net; % model from previous stage
    layers_ft = net_last.Layers;
end
[net, info] = trainNetwork(x_tra, y_tra, layers_ft, options);
```

|Condition|C1|C2|| | | | | |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
|Tra-RMSE /mm	|1.2510|1.3681| | | | | | |
|Tra-R2			|0.9945|0.9940| | | | | | |
|Tra-Bias /mm	|0.00|-0.04| | | | | | |
|Tra-95CI /mm	|[-2.45, 2.46]|[-2.72, 2.64]| | | | | | |
|Tes-RMSE /mm	|3.0608|1.4264| | | | | | |
|Tes-R2			|0.9649|0.9935| | | | | | |
|Tes-Bias /mm	|0.45|-0.04| | | | | | |
|Tes-95CI /mm	|[-5.48, 6.39]|[-2.83, 2.76]| | | | | | |

	- C1: previous model trained using the optimal architecture and paras and only left or right data
 	- C2: model trained using both left and right data and shuffled data windows 

|Condition|C3|C4|C5|C6|C7|
|:---|---:|---:|---:|---:|---:|
|Norm Test-RMSE /mm	|1.3555|1.7794|1.8207|1.8578|1.8334|
|Norm Test-R2		|0.9942|0.9900|0.9893|0.9890|0.9888|
|Norm Test-Bias /mm	|-0.05|-0.21|-0.19|-0.19|-0.16|
|Norm Test-95CI /mm	|[-2.70, 2.61]|[-3.67, 3.26]|[-3.74, 3.36]|[-3.81, 3.43]|[-3.74, 3.42]|
|Edge Test-RMSE /mm	|52.9486|9.3420|5.0955|4.7546|3.2379|
|Edge Test-R2		|0.2795|0.9776|0.9933|0.9947|0.9970|
|Edge Test-Bias /mm	|-9.24|1.02|-0.19|-0.78|-0.01|
|Edge Test-95CI /mm	|[-111.43, 92.95]|[-17.19, 19.22]|[-10.17, 9.79]|[-9.97, 8.41]|[-6.36, 6.33]|
|C2 Test-RMSE /mm	|1.4264|1.8070|1.8555|1.9269|1.9072|
|C2 Test-R2			|0.9935|0.9895|0.9889|0.9881|0.9883|
|C2 Test-Bias /mm	|-0.04|-0.19|-0.17|-0.16|-0.14|
|C2 Test-95CI /mm	|[-2.83, 2.76]|[-3.71, 3.33]|[-3.79, 3.45]|[-3.92, 3.60]|[-3.87, 3.58]|

|Condition|C8|C9|C10|
|:---|---:|---:|---:|
|Norm Test-RMSE /mm	|1.9115|1.8170|1.9770|
|Norm Test-R2		|0.9888|0.9888|0.9875|
|Norm Test-Bias /mm	|0.10|-0.03|0.13|
|Norm Test-95CI /mm	|[-3.64, 3.84]|[-3.59, 3.53]|[-3.74, 4.00]|
|Edge Test-RMSE /mm	|3.2565|2.6616|2.3984|
|Edge Test-R2		|0.9973|0.9978|0.9987|
|Edge Test-Bias /mm	|-0.42|-0.06|-0.19|
|Edge Test-95CI /mm	|[-6.75, 5.91]|[-5.27, 5.16]|[-4.87, 4.50]|
|C2 Test-RMSE /mm	|1.8803|1.9223|1.9625|
|C2 Test-R2			|0.9886|0.9881|0.9876|
|C2 Test-Bias /mm	|0.08|-0.03|0.13|
|C2 Test-95CI /mm	|[-3.60, 3.76]|[-3.80, 3.74]|[-3.71, 3.97]|
	- 'Test' indicates the data used to final test the model performance.
 	- 'C2 Test' indicates the test data used in C2 model trial. Note that the size of C2 Test data is much larger than the Norm or Edge Test data, since Norm Test size had to be considerably reduced to achieve the required edge-to-normal ratio as the edge proportion increased. C2 Test data were never used during the fine-tunning process, so they are perfect to evaluate the degradation of model performance after fine-tunning.
	- C3: model without fine tunning
 	- C4: based on the model from C3, fine-tunned with a dataset containing 5% edge cases
 	- C5: based on the model from C4, fine-tunned with a dataset containing 10% edge cases
 	- C6: based on the model from C5, fine-tunned with a dataset containing 15% edge cases
 	- C7: based on the model from C6, fine-tunned with a dataset containing 20% edge cases
   	- C8: based on the model from C7, fine-tunned with a dataset containing 25% edge cases
   	- C9: based on the model from C8, fine-tunned with a dataset containing 30% edge cases
   	- C10: based on the model from C9, fine-tunned with a dataset containing 35% edge cases
