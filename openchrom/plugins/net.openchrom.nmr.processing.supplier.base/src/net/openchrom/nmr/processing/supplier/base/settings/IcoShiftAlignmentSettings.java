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
package net.openchrom.nmr.processing.supplier.base.settings;

import java.util.List;

import net.openchrom.nmr.processing.supplier.base.settings.support.IcoShiftAlignmentGapFillingType;
import net.openchrom.nmr.processing.supplier.base.settings.support.IcoShiftAlignmentShiftCorrectionType;
import net.openchrom.nmr.processing.supplier.base.settings.support.IcoShiftAlignmentTargetCalculationSelection;
import net.openchrom.nmr.processing.supplier.base.settings.support.IcoShiftAlignmentType;
import net.openchrom.nmr.processing.supplier.base.settings.support.IcoShiftAlignmentUtilities.ChemicalShiftInterval;

public class IcoShiftAlignmentSettings {

	private IcoShiftAlignmentTargetCalculationSelection targetCalculationSelection = IcoShiftAlignmentTargetCalculationSelection.MEAN;
	private IcoShiftAlignmentShiftCorrectionType shiftCorrectionType = IcoShiftAlignmentShiftCorrectionType.FAST;
	private int shiftCorrectionTypeValue;
	private IcoShiftAlignmentGapFillingType gapFillingType = IcoShiftAlignmentGapFillingType.MARGIN;
	private IcoShiftAlignmentType alignmentType = IcoShiftAlignmentType.SINGLE_PEAK;
	private double singlePeakLowerBorder;
	private double singlePeakHigherBorder;
	private int numberOfIntervals;
	private int intervalLength;
	private List<ChemicalShiftInterval> userDefIntervalRegions;
	private boolean preliminaryCoShifting;

	public IcoShiftAlignmentTargetCalculationSelection getTargetCalculationSelection() {

		return targetCalculationSelection;
	}

	public void setTargetCalculationSelection(IcoShiftAlignmentTargetCalculationSelection targetCalculationSelection) {

		this.targetCalculationSelection = targetCalculationSelection;
	}

	public IcoShiftAlignmentShiftCorrectionType getShiftCorrectionType() {

		return shiftCorrectionType;
	}

	public void setShiftCorrectionType(IcoShiftAlignmentShiftCorrectionType shiftCorrectionType) {

		this.shiftCorrectionType = shiftCorrectionType;
	}

	public IcoShiftAlignmentGapFillingType getGapFillingType() {

		return gapFillingType;
	}

	public void setGapFillingType(IcoShiftAlignmentGapFillingType gapFillingType) {

		this.gapFillingType = gapFillingType;
	}

	public IcoShiftAlignmentType getAlignmentType() {

		return alignmentType;
	}

	public void setAlignmentType(IcoShiftAlignmentType alignmentType) {

		this.alignmentType = alignmentType;
	}

	public double getSinglePeakLowerBorder() {

		return singlePeakLowerBorder;
	}

	public void setSinglePeakLowerBorder(double singlePeakLowerBorder) {

		this.singlePeakLowerBorder = singlePeakLowerBorder;
	}

	public double getSinglePeakHigherBorder() {

		return singlePeakHigherBorder;
	}

	public void setSinglePeakHigherBorder(double singlePeakHigherBorder) {

		this.singlePeakHigherBorder = singlePeakHigherBorder;
	}

	public int getIntervalLength() {

		return intervalLength;
	}

	public void setIntervalLength(int intervalLength) {

		this.intervalLength = intervalLength;
	}

	public List<ChemicalShiftInterval> getUserDefIntervalRegions() {

		return userDefIntervalRegions;
	}

	public void setUserDefIntervalRegions(List<ChemicalShiftInterval> userDefIntervalRegions) {

		this.userDefIntervalRegions = userDefIntervalRegions;
	}

	public int getNumberOfIntervals() {

		return numberOfIntervals;
	}

	public void setNumberOfIntervals(int numberOfIntrvals) {

		this.numberOfIntervals = numberOfIntrvals;
	}

	public int getShiftCorrectionTypeValue() {

		return shiftCorrectionTypeValue;
	}

	public void setShiftCorrectionTypeValue(int shiftCorrectionTypeValue) {

		this.shiftCorrectionTypeValue = shiftCorrectionTypeValue;
	}

	public boolean isPreliminaryCoShifting() {

		return preliminaryCoShifting;
	}

	public void setPreliminaryCoShifting(boolean preliminaryCoShifting) {

		this.preliminaryCoShifting = preliminaryCoShifting;
	}
}