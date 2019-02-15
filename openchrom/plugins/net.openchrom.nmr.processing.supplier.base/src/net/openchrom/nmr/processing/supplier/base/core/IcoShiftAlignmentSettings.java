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

import java.util.List;

import net.openchrom.nmr.processing.supplier.base.core.IcoShiftAlignment.AlignmentType;
import net.openchrom.nmr.processing.supplier.base.core.IcoShiftAlignment.ChemicalShiftInterval;
import net.openchrom.nmr.processing.supplier.base.core.IcoShiftAlignment.GapFillingType;
import net.openchrom.nmr.processing.supplier.base.core.IcoShiftAlignment.PrelimiterCoShifting;
import net.openchrom.nmr.processing.supplier.base.core.IcoShiftAlignment.ShiftCorrectionType;
import net.openchrom.nmr.processing.supplier.base.core.IcoShiftAlignment.TargetCalculationSelection;

public class IcoShiftAlignmentSettings {

	private TargetCalculationSelection targetCalculationSelection = TargetCalculationSelection.MEAN;
	private ShiftCorrectionType shiftCorrectionType = ShiftCorrectionType.FAST;
	private int shiftCorrectionTypeValue;
	private GapFillingType gapFillingType = GapFillingType.ZERO;
	private AlignmentType alignmentType = AlignmentType.SINGLE_PEAK;
	private double singlePeakLowerBorder;
	private double singlePeakHigherBorder;
	private int numberOfIntervals;
	private double intervalLength;
	private List<ChemicalShiftInterval> userDefIntervalRegions;
	private PrelimiterCoShifting prelimiterCoShifting = PrelimiterCoShifting.NO;

	public TargetCalculationSelection getTargetCalculationSelection() {

		return targetCalculationSelection;
	}

	public void setTargetCalculationSelection(TargetCalculationSelection targetCalculationSelection) {

		this.targetCalculationSelection = targetCalculationSelection;
	}

	public ShiftCorrectionType getShiftCorrectionType() {

		return shiftCorrectionType;
	}

	public void setShiftCorrectionType(ShiftCorrectionType shiftCorrectionType) {

		this.shiftCorrectionType = shiftCorrectionType;
	}

	public GapFillingType getGapFillingType() {

		return gapFillingType;
	}

	public void setGapFillingType(GapFillingType gapFillingType) {

		this.gapFillingType = gapFillingType;
	}

	public AlignmentType getAlignmentType() {

		return alignmentType;
	}

	public void setAligmentType(AlignmentType intervalSelection) {

		this.alignmentType = intervalSelection;
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

	public double getIntervalLength() {

		return intervalLength;
	}

	public void setIntervalLength(double intervalLength) {

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

	public void setPrelimiterCoShifting(PrelimiterCoShifting prelimiterCoShifting) {

		this.prelimiterCoShifting = prelimiterCoShifting;
	}

	public PrelimiterCoShifting getPrelimiterCoShifting() {

		return prelimiterCoShifting;
	}
}
