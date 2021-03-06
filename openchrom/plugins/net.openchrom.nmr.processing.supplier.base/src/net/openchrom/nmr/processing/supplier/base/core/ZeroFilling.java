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

import org.apache.commons.math3.complex.Complex;
import org.eclipse.chemclipse.nmr.model.core.IScanNMR;
import org.eclipse.chemclipse.nmr.processor.core.AbstractScanProcessor;
import org.eclipse.chemclipse.nmr.processor.core.IScanProcessor;
import org.eclipse.chemclipse.nmr.processor.settings.IProcessorSettings;
import org.eclipse.chemclipse.processing.core.IProcessingInfo;
import org.eclipse.core.runtime.IProgressMonitor;

import net.openchrom.nmr.processing.supplier.base.settings.ZeroFillingSettings;

public class ZeroFilling extends AbstractScanProcessor implements IScanProcessor {

	@Override
	public IProcessingInfo process(IScanNMR scanNMR, IProcessorSettings processorSettings, IProgressMonitor monitor) {

		IProcessingInfo processingInfo = validate(scanNMR, processorSettings, ZeroFillingSettings.class);
		if(!processingInfo.hasErrorMessages()) {
			ZeroFillingSettings settings = (ZeroFillingSettings)processorSettings;
			Complex[] modifiedSignals = scanNMR.getModifiedSignals();
			/*
			 * Add your implementation here.
			 */
		}
		return processingInfo;
	}
}
