/*******************************************************************************
 * Copyright (c) 2018 Lablicate GmbH.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 * 
 * Contributors:
 * jan - initial API and implementation
 *******************************************************************************/
package net.openchrom.nmr.processing.supplier.base.core;

import java.util.NavigableSet;

import org.eclipse.chemclipse.nmr.model.core.IScanNMR;
import org.eclipse.chemclipse.nmr.model.core.ISignalFID;
import org.eclipse.chemclipse.nmr.processor.core.AbstractScanProcessor;
import org.eclipse.chemclipse.nmr.processor.core.IScanProcessor;
import org.eclipse.chemclipse.nmr.processor.settings.IProcessorSettings;
import org.eclipse.chemclipse.processing.core.IProcessingInfo;
import org.eclipse.core.runtime.IProgressMonitor;

public class ExponentialApodizationFunction extends AbstractScanProcessor implements IScanProcessor {

	@Override
	public IProcessingInfo process(IScanNMR scanNMR, IProcessorSettings processorSettings, IProgressMonitor monitor) {

		final IProcessingInfo processingInfo = validate(scanNMR, processorSettings);
		if(!processingInfo.hasErrorMessages()) {
			NavigableSet<ISignalFID> signals = scanNMR.getSignalsFID();
			double time[] = signals.stream().mapToDouble(ISignalFID::getTime).toArray();
			double[] data = exponentialApodizationFunction(time, scanNMR);
			int i = 0;
			for(ISignalFID signal : signals) {
				signal.setIntensity(data[i] * signal.getIntensity());
				i++;
			}
			processingInfo.setProcessingResult(scanNMR);
		}
		return processingInfo;
	}

	private double[] exponentialApodizationFunction(double[] timeScale, IScanNMR scanNMR) {

		double exponentialLineBroadeningFactor = 0;
		if(scanNMR.processingParametersContainsKey("exponentialLineBroadeningFactor")) {
			exponentialLineBroadeningFactor = scanNMR.getProcessingParameters("exponentialLineBroadeningFactor");
		}
		double[] exponentialLineBroadening = new double[timeScale.length];
		double exponentialLineBroadenigTerm;
		if(exponentialLineBroadeningFactor > 0) {
			for(int i = 0; i < timeScale.length; i++) { // Lbfunc=exp(-Timescale'*pi*NmrData.lb);
				exponentialLineBroadenigTerm = (-timeScale[i] * Math.PI * exponentialLineBroadeningFactor);
				exponentialLineBroadening[i] = Math.exp(exponentialLineBroadenigTerm);
			}
		} else {
			for(int i = 0; i < timeScale.length; i++) {
				exponentialLineBroadening[i] = (timeScale[i] * 0 + 1);
			}
		}
		return exponentialLineBroadening;
	}
}
