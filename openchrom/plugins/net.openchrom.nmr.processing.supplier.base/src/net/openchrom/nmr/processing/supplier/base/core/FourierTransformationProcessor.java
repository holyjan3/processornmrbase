/*******************************************************************************
 * Copyright (c) 2018 Lablicate GmbH.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 * Dr. Philip Wenig - initial API and implementation
 *******************************************************************************/
package net.openchrom.nmr.processing.supplier.base.core;

import java.util.Arrays;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.eclipse.chemclipse.nmr.model.core.IScanNMR;
import org.eclipse.chemclipse.nmr.model.core.ISignalFID;
import org.eclipse.chemclipse.nmr.model.support.ISignalExtractor;
import org.eclipse.chemclipse.nmr.model.support.SignalExtractor;
import org.eclipse.chemclipse.nmr.processor.core.AbstractScanProcessor;
import org.eclipse.chemclipse.nmr.processor.core.IScanProcessor;
import org.eclipse.chemclipse.nmr.processor.settings.IProcessorSettings;
import org.eclipse.chemclipse.processing.core.IProcessingInfo;
import org.eclipse.core.runtime.IProgressMonitor;

import net.openchrom.nmr.processing.supplier.base.settings.FourierTransformationSettings;

public class FourierTransformationProcessor extends AbstractScanProcessor implements IScanProcessor {

	@Override
	public IProcessingInfo process(final IScanNMR scanNMR, final IProcessorSettings processorSettings, final IProgressMonitor monitor) {

		IProcessingInfo processingInfo = validate(scanNMR, processorSettings);
		if(!processingInfo.hasErrorMessages()) {
			FourierTransformationSettings settings = (FourierTransformationSettings)processorSettings;
			ISignalExtractor signalExtractor = new SignalExtractor(scanNMR);
			Complex[] fourierTransformedData = transform(scanNMR, settings);
			double[] chemicalShift = generateChemicalShiftAxis(scanNMR);
			signalExtractor.createScans(fourierTransformedData, chemicalShift);
			processingInfo.setProcessingResult(scanNMR);
		}
		return processingInfo;
	}

	public static double[] generateChemicalShiftAxis(IScanNMR scanNMR) {

		double doubleSize = scanNMR.getProcessingParameters("numberOfFourierPoints");
		int deltaAxisPoints = (int)doubleSize;
		double[] chemicalShiftAxis = new double[(int)doubleSize];
		double minValueDeltaAxis = scanNMR.getProcessingParameters("firstDataPointOffset");
		double maxValueDeltaAxis = scanNMR.getProcessingParameters("sweepWidth") + scanNMR.getProcessingParameters("firstDataPointOffset");
		UtilityFunctions utilityFunction = new UtilityFunctions();
		chemicalShiftAxis = utilityFunction.generateLinearlySpacedVector(minValueDeltaAxis, maxValueDeltaAxis, deltaAxisPoints);
		return chemicalShiftAxis;
	}

	private Complex[] transform(IScanNMR scanNMR, FourierTransformationSettings processorSettings) {

		/*
		 * Raw Data
		 */
		double[] signals = scanNMR.getSignalsFID().stream().mapToDouble(ISignalFID::getIntensity).toArray();
		Complex[] freeInductionDecayShiftedWindowMultiplication = Arrays.stream(signals).boxed().map(Complex::new).toArray(Complex[]::new);
		// zero filling // Automatic zero filling if size != 2^n
		int zeroFillingFactor = 0;
		Complex[] freeInductionDecayShiftedWindowMultiplicationZeroFill = new Complex[freeInductionDecayShiftedWindowMultiplication.length];
		int checkPowerOfTwo = freeInductionDecayShiftedWindowMultiplication.length % 256;
		boolean automaticZeroFill = true;
		if(checkPowerOfTwo > 0) {
			for(int i = 10; i < 17; i++) {
				int automaticSize = (int)Math.pow(2, i);
				if(automaticSize > freeInductionDecayShiftedWindowMultiplication.length) {
					freeInductionDecayShiftedWindowMultiplicationZeroFill = new Complex[automaticSize];
					for(int j = 0; j < automaticSize; j++) {
						freeInductionDecayShiftedWindowMultiplicationZeroFill[j] = new Complex(0, 0);
					}
					int copySize = freeInductionDecayShiftedWindowMultiplication.length;
					System.arraycopy(freeInductionDecayShiftedWindowMultiplication, 0, freeInductionDecayShiftedWindowMultiplicationZeroFill, 0, copySize);
					int numberOfFourierPoints = automaticSize;
					scanNMR.putProcessingParameters("numberOfFourierPoints", (double)numberOfFourierPoints);
					automaticZeroFill = false;
					break;
				}
			}
		}
		if(zeroFillingFactor == 1) { // 16k
			int newDataSize = (int)Math.pow(2, 14);
			freeInductionDecayShiftedWindowMultiplicationZeroFill = new Complex[newDataSize];
			for(int i = 0; i < newDataSize; i++) {
				freeInductionDecayShiftedWindowMultiplicationZeroFill[i] = new Complex(0, 0);
			}
			int copySize = freeInductionDecayShiftedWindowMultiplication.length;
			System.arraycopy(freeInductionDecayShiftedWindowMultiplication, 0, freeInductionDecayShiftedWindowMultiplicationZeroFill, 0, copySize);
			double numberOfFourierPoints = newDataSize;
			scanNMR.putProcessingParameters("numberOfFourierPoints", numberOfFourierPoints);
		} else if(zeroFillingFactor == 2) { // 32k
			int newDataSize = (int)Math.pow(2, 15);
			freeInductionDecayShiftedWindowMultiplicationZeroFill = new Complex[newDataSize];
			for(int i = 0; i < newDataSize; i++) {
				freeInductionDecayShiftedWindowMultiplicationZeroFill[i] = new Complex(0, 0);
			}
			int copySize = freeInductionDecayShiftedWindowMultiplication.length;
			System.arraycopy(freeInductionDecayShiftedWindowMultiplication, 0, freeInductionDecayShiftedWindowMultiplicationZeroFill, 0, copySize);
			double numberOfFourierPoints = newDataSize;
			scanNMR.putProcessingParameters("numberOfFourierPoints", numberOfFourierPoints);
		} else if(zeroFillingFactor == 3) { // 64k
			int newDataSize = (int)Math.pow(2, 16);
			freeInductionDecayShiftedWindowMultiplicationZeroFill = new Complex[newDataSize];
			for(int i = 0; i < newDataSize; i++) {
				freeInductionDecayShiftedWindowMultiplicationZeroFill[i] = new Complex(0, 0);
			}
			int copySize = freeInductionDecayShiftedWindowMultiplication.length;
			System.arraycopy(freeInductionDecayShiftedWindowMultiplication, 0, freeInductionDecayShiftedWindowMultiplicationZeroFill, 0, copySize);
			double numberOfFourierPoints = newDataSize;
			scanNMR.putProcessingParameters("numberOfFourierPoints", numberOfFourierPoints);
		} else {
			// do nothing
			if(automaticZeroFill) {
				int dataSize = freeInductionDecayShiftedWindowMultiplication.length;
				freeInductionDecayShiftedWindowMultiplicationZeroFill = new Complex[dataSize];
				System.arraycopy(freeInductionDecayShiftedWindowMultiplication, 0, freeInductionDecayShiftedWindowMultiplicationZeroFill, 0, dataSize);
				double numberOfFourierPoints = dataSize;
				scanNMR.putProcessingParameters("numberOfFourierPoints", numberOfFourierPoints);
			}
		}
		// // generate x-axis (delta [ppm])
		// double[] deltaAxisPPM = generateChemicalShiftAxis(brukerParameterMap);
		// Fourier transform, shift and flip the data
		Complex[] nmrSpectrumProcessed = fourierTransformNmrData(freeInductionDecayShiftedWindowMultiplicationZeroFill);
		if(scanNMR.getProcessingParameters("ProcessedDataFlag").equals(1.0)) {
			// shift processed data once more
			leftShiftNMRComplexData(nmrSpectrumProcessed, nmrSpectrumProcessed.length / 2);
		}
		return nmrSpectrumProcessed;
	}

	@SuppressWarnings("unused")
	private void leftShiftNMRData(int[] dataArray, int pointsToShift) {

		pointsToShift = pointsToShift % dataArray.length;
		while(pointsToShift-- > 0) {
			int tempArray = dataArray[0];
			for(int i = 1; i < dataArray.length; i++) {
				dataArray[i - 1] = dataArray[i];
			}
			dataArray[dataArray.length - 1] = tempArray;
		}
	}

	private int[] rightShiftNMRData(int[] dataArray, int pointsToShift) {

		for(int i = 0; i < pointsToShift; i++) {
			int tempArray = dataArray[dataArray.length - 1];
			for(int g = dataArray.length - 2; g > -1; g--) {
				dataArray[g + 1] = dataArray[g];
			}
			dataArray[0] = tempArray;
		}
		return dataArray;
	}

	private void leftShiftNMRComplexData(Complex[] dataArray, int pointsToShift) {

		pointsToShift = pointsToShift % dataArray.length;
		while(pointsToShift-- > 0) {
			Complex tempArray = dataArray[0];
			for(int i = 1; i < dataArray.length; i++) {
				dataArray[i - 1] = dataArray[i];
			}
			dataArray[dataArray.length - 1] = tempArray;
		}
	}

	private Complex[] rightShiftNMRComplexData(Complex[] dataArray, int pointsToShift) {

		for(int i = 0; i < pointsToShift; i++) {
			Complex tempArray = dataArray[dataArray.length - 1];
			for(int g = dataArray.length - 2; g > -1; g--) {
				dataArray[g + 1] = dataArray[g];
			}
			dataArray[0] = tempArray;
		}
		return dataArray;
	}

	private Complex[] fourierTransformNmrData(Complex[] fid) {

		FastFourierTransformer fFourierTransformer = new FastFourierTransformer(DftNormalization.STANDARD);
		Complex[] nmrSpectrum = fFourierTransformer.transform(fid, TransformType.FORWARD);
		Complex[] nmrSpectrumProcessed = new Complex[nmrSpectrum.length];
		System.arraycopy(nmrSpectrum, 0, nmrSpectrumProcessed, 0, nmrSpectrum.length); // NmrData.SPECTRA
		rightShiftNMRComplexData(nmrSpectrumProcessed, nmrSpectrumProcessed.length / 2);
		ArrayUtils.reverse(nmrSpectrumProcessed);
		return nmrSpectrumProcessed;
	}
}
