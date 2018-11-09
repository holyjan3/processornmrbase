/*******************************************************************************
 * Copyright (c) 2018 Lablicate GmbH.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 * 
 * Contributors:
 * Jan Holy - initial API and implementation
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

public class GaussianApodizationFunction extends AbstractScanProcessor implements IScanProcessor {

	@Override
	public IProcessingInfo process(IScanNMR scanNMR, IProcessorSettings processorSettings, IProgressMonitor monitor) {

		final IProcessingInfo processingInfo = validate(scanNMR, processorSettings);
		if(!processingInfo.hasErrorMessages()) {
			NavigableSet<ISignalFID> signals = scanNMR.getSignalsFID();
			double time[] = signals.stream().mapToDouble(ISignalFID::getTime).toArray();
			double[] data = gaussianApodizationFunction(time, scanNMR);
			int i = 0;
			for(ISignalFID signal : signals) {
				signal.setIntensity(data[i] * signal.getIntensity());
				i++;
			}
			processingInfo.setProcessingResult(scanNMR);
		}
		return processingInfo;
	}

	private double[] gaussianApodizationFunction(double[] timeScale, IScanNMR scanNMR) {

		double gaussianLineBroadeningFactor = 0;
		if(scanNMR.processingParametersContainsKey("gaussianLineBroadeningFactor")) {
			gaussianLineBroadeningFactor = scanNMR.getProcessingParameters("gaussianLineBroadeningFactor");
		}
		double[] gaussianLineBroadening = new double[timeScale.length];
		double gaussianLineBroadenigTermA;
		double gaussianLineBroadenigTermB;
		if(gaussianLineBroadeningFactor > 0) {
			// gf=2*sqrt(log(2))/(pi*NmrData.gw);
			// Gwfunc=exp(-(Timescale'/gf).^2);
			gaussianLineBroadenigTermA = (Math.PI * gaussianLineBroadeningFactor);
			double gaussFactor = 2 * Math.sqrt(Math.log(2)) / gaussianLineBroadenigTermA;
			for(int i = 0; i < timeScale.length; i++) {
				gaussianLineBroadenigTermB = -(timeScale[i] / gaussFactor);
				gaussianLineBroadening[i] = Math.exp(Math.pow(gaussianLineBroadenigTermB, 2));
			}
		} else {
			for(int i = 0; i < timeScale.length; i++) {
				gaussianLineBroadening[i] = (timeScale[i] * 0 + 1);
			}
		}
		return gaussianLineBroadening;
	}
}
