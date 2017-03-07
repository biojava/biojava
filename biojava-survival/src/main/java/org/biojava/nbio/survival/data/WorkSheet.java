/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.survival.data;

import java.io.*;
import java.util.*;

/**
 * Need to handle very large spreadsheets of expression data so keep memory
 * footprint low
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class WorkSheet {

	private LinkedHashMap<String, HeaderInfo> columnLookup = new LinkedHashMap<String, HeaderInfo>();
	private LinkedHashMap<String, HeaderInfo> rowLookup = new LinkedHashMap<String, HeaderInfo>();
	private CompactCharSequence[][] data = new CompactCharSequence[1][1];
	HashMap<String, String> dataGrid = new HashMap<String, String>();
	private String indexColumnName = "";

	/**
	 *
	 */
	public WorkSheet() {
	}

	/**
	 *
	 * @param rows
	 * @param columns
	 * @throws Exception
	 */
	public WorkSheet(Collection<String> rows, Collection<String> columns) throws Exception {
		//    rowsList = new ArrayList<String>(rows);
		int i = 1;
		for (String row : rows) {
			if (rowLookup.containsKey(row)) {
				throw new Exception("Duplicate row " + row);
			}
			rowLookup.put(row, new HeaderInfo(i));
			i++;
		}
		i = 1;
		for (String col : columns) {
			if (columnLookup.containsKey(col)) {
				throw new Exception("Duplicate row " + col);
			}
			columnLookup.put(col, new HeaderInfo(i));
			i++;
		}



		//  columnsList.trimToSize();
		//  rowsList.trimToSize();
		data = new CompactCharSequence[rowLookup.size() + 1][columnLookup.size() + 1];
	}

	/**
	 *
	 * @param values
	 */
	public WorkSheet(String[][] values) {
		//    System.out.println("In worksheet init " + Runtime.getRuntime().totalMemory());
		String[] columns = new String[values[0].length];
		for (int i = 0; i < columns.length; i++) {
			columns[i] = new String(values[0][i].getBytes());
		}
		for (int i = 1; i < columns.length; i++) {
			columnLookup.put(columns[i], new HeaderInfo(i));
		}


		for (int i = 1; i < values.length; i++) {
			String row = new String(values[i][0].getBytes());
			rowLookup.put(row, new HeaderInfo(i));
		}

		data = new CompactCharSequence[values.length][values[0].length];
		for (int row = 0; row < values.length; row++) {
			for (int col = 0; col < values[0].length; col++) {
				String value = values[row][col];
				data[row][col] = new CompactCharSequence(value);
				values[row][col] = null;
			}
			System.out.println("Row " + row + " " + Runtime.getRuntime().totalMemory());

		}
		values = null;
		System.gc();
		//data = values;

	}

	/**
	 * See if we can free up memory
	 */
	public void clear() {
		columnLookup.clear();
		rowLookup.clear();
		data = null;
		dataGrid.clear();
		doubleValues.clear();
		System.gc();
	}

	@Override
	public String toString() {
		return super.toString(); //To change body of generated methods, choose Tools | Templates.
	}

	/**
	 * Split a worksheet randomly. Used for creating a discovery/validation data
	 * set The first file name will matched the percentage and the second file
	 * the remainder
	 *
	 * @param percentage
	 * @param fileName1
	 * @param fileName2
	 * @throws Exception
	 */
	public void randomlyDivideSave(double percentage, String fileName1, String fileName2) throws Exception {
		ArrayList<String> rows = this.getDataRows();
		Collections.shuffle(rows);
		int portion = (int) (rows.size() * percentage);
		for (int i = 0; i < portion; i++) {
			this.hideRow(rows.get(i), true);
		}
		this.saveTXT(fileName2);
		for (int i = 0; i < portion; i++) {
			this.hideRow(rows.get(i), false);
		}
		for (int i = portion; i < rows.size(); i++) {
			this.hideRow(rows.get(i), true);
		}
		this.saveTXT(fileName1);
		for (int i = portion; i < rows.size(); i++) {
			this.hideRow(rows.get(i), false);
		}

	}

	/**
	 * Create a copy of a worksheet. If shuffling of columns or row for testing
	 * a way to duplicate original worksheet
	 *
	 * @param copyWorkSheet
	 * @param rows
	 * @return
	 * @throws Exception
	 */
	static public WorkSheet getCopyWorkSheetSelectedRows(WorkSheet copyWorkSheet, ArrayList<String> rows) throws Exception {

		ArrayList<String> columns = copyWorkSheet.getColumns();


		WorkSheet workSheet = new WorkSheet(rows, columns);
		for (String row : rows) {
			for (String col : columns) {
				workSheet.addCell(row, col, copyWorkSheet.getCell(row, col));
			}
		}
		workSheet.setMetaDataColumns(copyWorkSheet.getMetaDataColumns());
		workSheet.setMetaDataRows(copyWorkSheet.getMetaDataRows());
		return workSheet;

	}

	/**
	 * Create a copy of a worksheet. If shuffling of columns or row for testing
	 * a way to duplicate original worksheet
	 *
	 * @param copyWorkSheet
	 * @return
	 * @throws Exception
	 */
	static public WorkSheet getCopyWorkSheet(WorkSheet copyWorkSheet) throws Exception {
		ArrayList<String> rows = copyWorkSheet.getRows();
		ArrayList<String> columns = copyWorkSheet.getColumns();


		WorkSheet workSheet = new WorkSheet(rows, columns);
		for (String row : rows) {
			for (String col : columns) {
				workSheet.addCell(row, col, copyWorkSheet.getCell(row, col));
			}
		}
		workSheet.setMetaDataColumns(copyWorkSheet.getMetaDataColumns());
		workSheet.setMetaDataRows(copyWorkSheet.getMetaDataRows());
		return workSheet;

	}

	/**
	 *
	 * @param values
	 */
	public WorkSheet(CompactCharSequence[][] values) {
		//     System.out.println("In worksheet init " + Runtime.getRuntime().totalMemory());
		String[] columns = new String[values[0].length];
		for (int i = 0; i < columns.length; i++) {
			columns[i] = values[0][i].toString();
		}
		this.setIndexColumnName(columns[0]);
		for (int i = 1; i < columns.length; i++) {
			columnLookup.put(columns[i], new HeaderInfo(i));
		}


		for (int i = 1; i < values.length; i++) {
			String row = values[i][0].toString();
			rowLookup.put(row, new HeaderInfo(i));
		}

		data = values;
	}
	private LinkedHashMap<String, String> metaDataColumnsHashMap = new LinkedHashMap<String, String>();

	/**
	 *
	 * @return
	 */
	public ArrayList<String> getMetaDataColumns() {
		ArrayList<String> metaColumns = new ArrayList<String>();
		for (String key : metaDataColumnsHashMap.keySet()) {
			HeaderInfo hi = columnLookup.get(key);
			if (!hi.isHide()) {
				metaColumns.add(key);
			}
		}
		return metaColumns;
	}

	/**
	 *
	 * @return
	 */
	public ArrayList<String> getMetaDataRows() {
		ArrayList<String> metaRows = new ArrayList<String>();
		for (String key : metaDataRowsHashMap.keySet()) {
			HeaderInfo hi = rowLookup.get(key);
			if (!hi.isHide()) {
				metaRows.add(key);
			}
		}
		return metaRows;
	}

	/**
	 *
	 * @return
	 */
	public ArrayList<String> getDataColumns() {
		ArrayList<String> dataColumns = new ArrayList<String>();
		ArrayList<String> columns = this.getColumns();
		for (String column : columns) {
			if (!metaDataColumnsHashMap.containsKey(column)) {
				dataColumns.add(column);
			}
		}
		return dataColumns;
	}

	/**
	 * Randomly shuffle the columns and rows. Should be constrained to the same
	 * data type if not probably doesn't make any sense.
	 *
	 * @param columns
	 * @param rows
	 * @throws Exception
	 */
	public void shuffleColumnsAndThenRows(ArrayList<String> columns, ArrayList<String> rows) throws Exception {
		doubleValues.clear();

		for (String column : columns) { //shuffle all values in the column
			ArrayList<Integer> rowIndex = new ArrayList<Integer>();
			for (int i = 0; i < rows.size(); i++) {
				rowIndex.add(i);
			}
			Collections.shuffle(rowIndex);
			for (int i = 0; i < rows.size(); i++) {
				String row = rows.get(i);
				int randomIndex = rowIndex.get(i);
				String destinationRow = rows.get(randomIndex);


				String temp = this.getCell(destinationRow, column);
				String value = this.getCell(row, column);
				this.addCell(destinationRow, column, value);
				this.addCell(row, column, temp);
			}
		}

		for (String row : rows) {
			ArrayList<Integer> columnIndex = new ArrayList<Integer>();
			for (int i = 0; i < columns.size(); i++) {
				columnIndex.add(i);
			}
			Collections.shuffle(columnIndex);
			for (int i = 0; i < columns.size(); i++) {
				String column = columns.get(i);

				int randomIndex = columnIndex.get(i);
				String destinationCol = columns.get(randomIndex);


				String temp = this.getCell(row, destinationCol);
				String value = this.getCell(row, column);
				this.addCell(row, destinationCol, value);
				this.addCell(row, column, temp);
			}
		}


	}

	/**
	 * Need to shuffle column values to allow for randomized testing. The
	 * columns in the list will be shuffled together
	 *
	 * @param columns
	 * @throws Exception
	 */
	public void shuffleColumnValues(ArrayList<String> columns) throws Exception {
		doubleValues.clear();
		ArrayList<String> rows = this.getDataRows();
		for (String column : columns) { //shuffle all values in the column
			ArrayList<Integer> rowIndex = new ArrayList<Integer>();
			for (int i = 0; i < rows.size(); i++) {
				rowIndex.add(i);
			}
			Collections.shuffle(rowIndex);
			for (int i = 0; i < rows.size(); i++) {
				String row = rows.get(i);
				int randomIndex = rowIndex.get(i);
				String destinationRow = rows.get(randomIndex);


				String temp = this.getCell(destinationRow, column);
				String value = this.getCell(row, column);
				this.addCell(destinationRow, column, value);
				this.addCell(row, column, temp);
			}
		}

	}

	/**
	 * Need to shuffle rows values to allow for randomized testing. The rows in
	 * the list will be shuffled together
	 *
	 * @param rows
	 * @throws Exception
	 */
	public void shuffleRowValues(ArrayList<String> rows) throws Exception {
		doubleValues.clear();
		ArrayList<String> columns = this.getColumns();
		for (String row : rows) {
			ArrayList<Integer> columnIndex = new ArrayList<Integer>();
			for (int i = 0; i < columns.size(); i++) {
				columnIndex.add(i);
			}
			Collections.shuffle(columnIndex);

			for (int i = 0; i < columns.size(); i++) {
				String column = columns.get(i);
				int randomIndex = columnIndex.get(i);
				String destinationCol = columns.get(randomIndex);

				String temp = this.getCell(row, destinationCol);
				String value = this.getCell(row, column);
				this.addCell(row, destinationCol, value);
				this.addCell(row, column, temp);
			}
		}

	}

	/**
	 *
	 * @param value
	 */
	public void hideMetaDataColumns(boolean value) {
		ArrayList<String> metadataColumns = this.getMetaDataColumns();
		for (String column : metadataColumns) {
			this.hideColumn(column, value);
		}
	}

	/**
	 *
	 * @param value
	 */
	public void hideMetaDataRows(boolean value) {
		ArrayList<String> metadataRows = this.getMetaDataRows();
		for (String row : metadataRows) {
			this.hideRow(row, value);
		}
	}

	/**
	 *
	 */
	public void setMetaDataRowsAfterRow() {
		this.setMetaDataRowsAfterRow("META_DATA");
	}

	/**
	 *
	 */
	public void setMetaDataColumnsAfterColumn() {
		this.setMetaDataColumnsAfterColumn("META_DATA");
	}

	/**
	 *
	 * @param row
	 */
	public void setMetaDataRowsAfterRow(String row) {
		ArrayList<String> rows = this.getRows();
		boolean metarow = false;
		for (String r : rows) {
			if (r.equals(row) && !metarow) {
				metarow = true;
			}
			if (metarow) {
				this.markMetaDataRow(r);
			}
		}
	}

	/**
	 *
	 * @param column
	 */
	public void setMetaDataColumnsAfterColumn(String column) {
		ArrayList<String> cols = this.getColumns();
		boolean metacolumns = false;
		for (String col : cols) {
			if (col.equals(column) && !metacolumns) {
				metacolumns = true;
			}
			if (metacolumns) {
				this.markMetaDataColumn(col);
			}
		}


	}

	/**
	 * Clears existing meta data columns and sets new ones
	 *
	 * @param metaDataColumns
	 */
	public void setMetaDataColumns(ArrayList<String> metaDataColumns) {
		metaDataColumnsHashMap.clear();
		markMetaDataColumns(metaDataColumns);
	}

	/**
	 * marks columns as containing meta data
	 *
	 * @param metaDataColumns
	 */
	public void markMetaDataColumns(ArrayList<String> metaDataColumns) {
		for (String column : metaDataColumns) {
			metaDataColumnsHashMap.put(column, column);
		}
	}

	/**
	 *
	 * @param column
	 */
	public void markMetaDataColumn(String column) {
		metaDataColumnsHashMap.put(column, column);
	}

	/**
	 *
	 * @param column
	 * @return
	 */
	public boolean isMetaDataColumn(String column) {
		if (metaDataColumnsHashMap.get(column) == null) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 *
	 * @param row
	 * @return
	 */
	public boolean isMetaDataRow(String row) {
		if (metaDataRowsHashMap.get(row) == null) {
			return false;
		} else {
			return true;
		}
	}
	private LinkedHashMap<String, String> metaDataRowsHashMap = new LinkedHashMap<String, String>();

	/**
	 *
	 * @param row
	 */
	public void markMetaDataRow(String row) {
		metaDataRowsHashMap.put(row, row);
	}

	/**
	 *
	 * @param metaDataRows
	 */
	public void setMetaDataRows(ArrayList<String> metaDataRows) {
		metaDataRowsHashMap.clear();
		for (String row : metaDataRows) {
			metaDataRowsHashMap.put(row, row);
		}
	}

	/**
	 *
	 * @throws Exception
	 */
	public void hideEmptyRows() throws Exception {
		ArrayList<String> rows = this.getDataRows();
		ArrayList<String> columns = this.getDataColumns();
		for (String row : rows) {
			boolean emptyRow = true;
			for (String column : columns) {
				String value = this.getCell(row, column).trim();
				if (value.length() > 0) {
					emptyRow = false;
					break;
				}
			}
			if (emptyRow) {
				hideRow(row, true);
			}
		}

	}

	/**
	 *
	 * @throws Exception
	 */
	public void hideEmptyColumns() throws Exception {
		ArrayList<String> rows = this.getDataRows();
		ArrayList<String> columns = this.getDataColumns();
		for (String column : columns) {
			boolean emptyColumn = true;
			for (String row : rows) {
				String value = this.getCell(row, column).trim();
				if (value.length() > 0) {
					emptyColumn = false;
					break;
				}
			}
			if (emptyColumn) {
				hideColumn(column, true);
			}
		}

	}

	/**
	 *
	 * @param row
	 * @param hide
	 */
	public void hideRow(String row, boolean hide) {
		HeaderInfo rowInfo = rowLookup.get(row);
		rowInfo.setHide(hide);
	}

	/**
	 *
	 * @param column
	 * @param hide
	 */
	public void hideColumn(String column, boolean hide) {
		HeaderInfo colInfo = columnLookup.get(column);
		colInfo.setHide(hide);
	}

	/**
	 * Change values in a column where 0 = something and 1 = something different
	 *
	 * @param column
	 * @param values
	 * @throws Exception
	 */
	public void replaceColumnValues(String column, HashMap<String, String> values) throws Exception {
		for (String row : rowLookup.keySet()) {
			String oldValue = this.getCell(row, column);
			String newValue = values.get(oldValue);
			this.addCell(row, column, newValue);
		}

	}

	/**
	 * Apply filter to a column to change values from say numberic to nominal
	 * based on some range
	 *
	 * @param column
	 * @param changeValue
	 * @throws Exception
	 */
	public void applyColumnFilter(String column, ChangeValue changeValue) throws Exception {
		for (String row : rowLookup.keySet()) {
			String oldValue = this.getCell(row, column);
			String newValue = changeValue.change(oldValue);
			this.addCell(row, column, newValue);
		}
	}

	/**
	 *
	 * @param column
	 * @param defaultValue
	 */
	public void addColumn(String column, String defaultValue) {
		ArrayList<String> columns = new ArrayList<String>();
		columns.add(column);
		addColumns(columns, defaultValue);
	}

	/**
	 * Add columns to worksheet and set default value
	 *
	 * @param columns
	 * @param defaultValue
	 */
	public void addColumns(ArrayList<String> columns, String defaultValue) {
		CompactCharSequence dv = new CompactCharSequence(defaultValue);
		for (int i = 0; i < data.length; i++) {
			CompactCharSequence[] row = data[i];
			int oldrowlength = data[i].length;
			data[i] = (CompactCharSequence[]) resizeArray(row, oldrowlength + columns.size());
			for (int c = 0; c < columns.size(); c++) {
				data[i][oldrowlength + c] = dv;
			}
			if (i == 0) {
				for (int c = 0; c < columns.size(); c++) {
					String column = columns.get(c);
					data[0][oldrowlength + c] = new CompactCharSequence(column);
					columnLookup.put(column, new HeaderInfo(oldrowlength + c));
				}
			}
		}
		//   columnLookup.get("ZNF30");

		//     int startIndex = columnLookup.size() + 1;
		//     for (String column : columns) {
		//        if(column.equals("ttr")){
		//            int dummy = 1;
		//        }
		//        columnLookup.put(column, new HeaderInfo(startIndex));
		//        startIndex++;
		//    }


	}

	/**
	 *
	 * @param row
	 * @param defaultValue
	 */
	public void addRow(String row, String defaultValue) {
		ArrayList<String> rows = new ArrayList<String>();
		rows.add(row);
		addRows(rows, defaultValue);
	}

	/**
	 * Add rows to the worksheet and fill in default value
	 *
	 * @param rows
	 * @param defaultValue
	 */
	public void addRows(ArrayList<String> rows, String defaultValue) {
		CompactCharSequence dv = new CompactCharSequence(defaultValue);
		int oldlength = data.length;
		int numColumns = 0;
		if (data.length > 0 && data[0] != null) {
			numColumns = data[0].length;
		}
		data = (CompactCharSequence[][]) resizeArray(data, data.length + rows.size());
		for (int r = 0; r < rows.size(); r++) {
			data[oldlength + r] = new CompactCharSequence[numColumns];
			for (int c = 0; c < numColumns; c++) {
				data[oldlength + r][c] = dv;
			}
			data[oldlength + r][0] = new CompactCharSequence(rows.get(r));
			rowLookup.put(rows.get(r), new HeaderInfo(r + oldlength));
		}
	}

	/**
	 * Reallocates an array with a new size, and copies the contents of the old
	 * array to the new array.
	 *
	 * @param oldArray the old array, to be reallocated.
	 * @param newSize the new array size.
	 * @return A new array with the same contents.
	 */
	private static Object resizeArray(Object oldArray, int newSize) {
		int oldSize = java.lang.reflect.Array.getLength(oldArray);
		Class<?> elementType = oldArray.getClass().getComponentType();
		Object newArray = java.lang.reflect.Array.newInstance(
				elementType, newSize);
		int preserveLength = Math.min(oldSize, newSize);
		if (preserveLength > 0) {
			System.arraycopy(oldArray, 0, newArray, 0, preserveLength);
		}
		return newArray;
	}

	/**
	 * Add data to a cell
	 *
	 * @param row
	 * @param col
	 * @param value
	 * @throws Exception
	 */
	public void addCell(String row, String col, String value) throws Exception {
		HeaderInfo rowIndex = rowLookup.get(row);
		HeaderInfo colIndex = columnLookup.get(col);
		if (rowIndex == null) {
			throw new Exception("Row " + row + " not found in worksheet");
		}
		if (colIndex == null) {
			throw new Exception("Column " + col + " not found in worksheet");
		}


		data[rowIndex.getIndex()][colIndex.getIndex()] = new CompactCharSequence(value);
	}

	/**
	 *
	 * @param row
	 * @return
	 */
	public boolean isValidRow(String row) {
		HeaderInfo rowIndex = rowLookup.get(row);
		if (rowIndex == null) {
			for (String rowtable : rowLookup.keySet()) {
				if (row.equalsIgnoreCase(rowtable)) {

					return true;
				}
			}
			return false;
		} else {
			return true;
		}
	}

	/**
	 *
	 * @param col
	 * @return
	 */
	public boolean isValidColumn(String col) {
		HeaderInfo colIndex = columnLookup.get(col);
		if (colIndex == null) {
			for (String coltable : columnLookup.keySet()) {
				if (col.equalsIgnoreCase(coltable)) {

					return true;
				}
			}

			return false;


		} else {
			return true;
		}
	}
	//When we do gene signatures we ask for the same data value often. This method took up 50% of the time.
	HashMap<String, Double> doubleValues = new HashMap<String, Double>();
	boolean cacheDoubleValues = false;

	/**
	 *
	 * @param value
	 */
	public void setCacheDoubleValues(boolean value) {
		cacheDoubleValues = value;
	}

	/**
	 *
	 * @param row
	 * @param col
	 * @return
	 * @throws Exception
	 */
	public Double getCellDouble(String row, String col) throws Exception {
		if (cacheDoubleValues) {
			String key = row + ":" + col;

			Double v = doubleValues.get(key);
			if (v != null) {
				return v;
			}
			String value = getCell(row, col);

			try {
				v = Double.parseDouble(value);
			} catch (Exception e) {
			}
			doubleValues.put(key, v);
			return v;
		} else {
			Double v = null;
			String value = getCell(row, col);
			try {
				v = Double.parseDouble(value);
			} catch (Exception e) {
			}
			return v;
		}

	}

	/**
	 * Get cell value
	 *
	 * @param row
	 * @param col
	 * @return
	 * @throws Exception
	 */
	public String getCell(String row, String col) throws Exception {
		if (col.equals(this.getIndexColumnName())) {
			return row;
		}
		HeaderInfo rowIndex = rowLookup.get(row);
		HeaderInfo colIndex = columnLookup.get(col);

		if (rowIndex == null) {
			//allow for case insentive search
			for (String rowtable : rowLookup.keySet()) {
				if (row.equalsIgnoreCase(rowtable)) {
					rowIndex = rowLookup.get(rowtable);
					break;
				}
			}
			if (rowIndex == null) {
				throw new Exception("Row " + row + " not found in worksheet");
			}
		}
		if (colIndex == null) {
			//allow for case insentive search
			for (String coltable : columnLookup.keySet()) {
				if (col.equalsIgnoreCase(coltable)) {
					colIndex = columnLookup.get(coltable);
					break;
				}
			}
			if (colIndex == null) {
				throw new Exception("Column " + col + " not found in worksheet");
			}
		}

		CompactCharSequence ccs = data[rowIndex.getIndex()][colIndex.getIndex()];
		if (ccs != null) {
			return ccs.toString();
		} else {
			return "";
		}

		// return .toString();
	}

	/**
	 *
	 * @param changeValue
	 */
	public void changeRowHeader(ChangeValue changeValue) {
		ArrayList<String> rows = new ArrayList<String>(rowLookup.keySet());
		for (String row : rows) {
			String newRow = changeValue.change(row);
			HeaderInfo value = rowLookup.get(row);
			rowLookup.remove(row);
			rowLookup.put(newRow, value);
		}
	}

	/**
	 *
	 * @param changeValue
	 */
	public void changeColumnHeader(ChangeValue changeValue) {
		ArrayList<String> columns = new ArrayList<String>(columnLookup.keySet());
		for (String col : columns) {
			String newCol = changeValue.change(col);
			HeaderInfo value = columnLookup.get(col);
			columnLookup.remove(col);
			columnLookup.put(newCol, value);
		}
	}

	/**
	 *
	 * @param row
	 * @param newRow
	 * @throws Exception
	 */
	public void changeRowHeader(String row, String newRow) throws Exception {
		HeaderInfo value = rowLookup.get(row);
		if (value == null) {
			throw new Exception("Row not found " + row);
		}
		rowLookup.remove(row);
		rowLookup.put(newRow, value);
		if (this.isMetaDataRow(row)) {
			metaDataRowsHashMap.remove(row);
			metaDataRowsHashMap.put(newRow, newRow);
		}
	}

	/**
	 * Change the columns in the HashMap Key to the name of the value
	 *
	 * @param newColumnValues
	 * @throws Exception
	 */
	public void changeColumnsHeaders(LinkedHashMap<String, String> newColumnValues) throws Exception {
		for (String oldColumn : newColumnValues.keySet()) {
			String newColumn = newColumnValues.get(oldColumn);
			changeColumnHeader(oldColumn, newColumn);
		}

	}

	/**
	 *
	 * @param col
	 * @param newCol
	 * @throws Exception
	 */
	public void changeColumnHeader(String col, String newCol) throws Exception {
		HeaderInfo value = columnLookup.get(col);
		if (value == null) {
			throw new Exception("Column not found " + col);
		}
		columnLookup.remove(col);
		columnLookup.put(newCol, value);
		if (this.isMetaDataColumn(col)) {
			metaDataColumnsHashMap.remove(col);
			metaDataColumnsHashMap.put(newCol, newCol);
		}

	}

	/**
	 *
	 * @param column
	 * @return
	 * @throws Exception
	 */
	public Integer getColumnIndex(String column) throws Exception {
		HeaderInfo headerInfo = columnLookup.get(column);
		if (headerInfo == null) {
			throw new Exception("Column " + column + " not found");
		}
		return headerInfo.getIndex();
	}

	/**
	 *
	 * @param row
	 * @return
	 * @throws Exception
	 */
	public Integer getRowIndex(String row) throws Exception {
		HeaderInfo headerInfo = rowLookup.get(row);
		if (headerInfo == null) {
			throw new Exception("Row " + row + " not found");
		}
		return headerInfo.getIndex();
	}

	/**
	 *
	 * @param number
	 * @return
	 */
	public ArrayList<String> getRandomDataColumns(int number) {
		ArrayList<String> columns = getDataColumns();
		return getRandomDataColumns(number, columns);
	}

	/**
	 *
	 * @param number
	 * @param columns
	 * @return
	 */
	public ArrayList<String> getRandomDataColumns(int number, ArrayList<String> columns) {
		ArrayList<String> randomColumns = new ArrayList<String>();
		HashMap<String, String> picked = new HashMap<String, String>();
		while (picked.size() < number) {
			double v = Math.random();
			int index = (int) (v * columns.size());
			if (picked.containsKey(String.valueOf(index))) {
				continue;
			}
			picked.put(String.valueOf(index), String.valueOf(index));
			randomColumns.add(columns.get(index));
		}
		return randomColumns;

	}

	/**
	 * Get the list of column names including those that may be hidden
	 *
	 * @return
	 */
	public ArrayList<String> getAllColumns() {
		ArrayList<String> columns = new ArrayList<String>();
		for (String col : columnLookup.keySet()) {
			columns.add(col);
		}
		return columns;
	}

	/**
	 * Get the list of column names. Does not include hidden columns
	 *
	 * @return
	 */
	public ArrayList<String> getColumns() {
		ArrayList<String> columns = new ArrayList<String>();
		for (String col : columnLookup.keySet()) {
			HeaderInfo hi = columnLookup.get(col);
			if (!hi.isHide()) {
				columns.add(col);
			}
		}
		return columns;
	}

	/**
	 * Get back a list of unique values in the column
	 *
	 * @param column
	 * @return
	 * @throws Exception
	 */
	public ArrayList<String> getDiscreteColumnValues(String column) throws Exception {
		HashMap<String, String> hashMapValues = new HashMap<String, String>();
		ArrayList<String> values = new ArrayList<String>();
		ArrayList<String> rows = getDataRows();
		for (String row : rows) {
			String value = getCell(row, column);
			if (!hashMapValues.containsKey(value)) {
				hashMapValues.put(value, value);
				values.add(value);
			}
		}
		return values;
	}

	/**
	 * Get back a list of unique values in the row
	 *
	 * @param row
	 * @return
	 * @throws Exception
	 */
	public ArrayList<String> getDiscreteRowValues(String row) throws Exception {
		HashMap<String, String> hashMapValues = new HashMap<String, String>();
		ArrayList<String> values = new ArrayList<String>();
		for (String column : getColumns()) {
			String value = getCell(row, column);
			if (!hashMapValues.containsKey(value)) {
				hashMapValues.put(value, value);
				values.add(value);
			}
		}
		return values;
	}

	/**
	 * Get all rows including those that may be hidden
	 *
	 * @return
	 */
	public ArrayList<String> getAllRows() {
		ArrayList<String> rows = new ArrayList<String>();
		for (String row : rowLookup.keySet()) {
			rows.add(row);
		}
		return rows;

	}

	/**
	 * Get the list of row names. Will exclude hidden values
	 *
	 * @return
	 */
	public ArrayList<String> getRows() {
		ArrayList<String> rows = new ArrayList<String>();
		for (String row : rowLookup.keySet()) {
			HeaderInfo hi = rowLookup.get(row);
			if (!hi.isHide()) {
				rows.add(row);
			}
		}
		return rows;
	}

	/**
	 * Get the list of row names
	 *
	 * @return
	 */
	public ArrayList<String> getDataRows() {
		ArrayList<String> rows = new ArrayList<String>();
		for (String row : rowLookup.keySet()) {
			if (this.isMetaDataRow(row)) {
				continue;
			}
			HeaderInfo hi = rowLookup.get(row);
			if (!hi.isHide()) {
				rows.add(row);
			}
		}
		return rows;
	}

	/**
	 * Get the log scale of this worksheet where a zero value will be set to .1
	 * as Log(0) is undefined
	 *
	 * @param base
	 * @return
	 * @throws Exception
	 */
	public WorkSheet getLogScale(double base) throws Exception {
		return getLogScale(base, .1);

	}

	/**
	 * Get the log scale of this worksheet
	 *
	 * @param base
	 * @return
	 * @throws Exception
	 */
	public WorkSheet getLogScale(double base, double zeroValue) throws Exception {

		WorkSheet workSheet = new WorkSheet(getRows(), getColumns());
		workSheet.setIndexColumnName(this.getIndexColumnName());
		ArrayList<String> rows = getRows();
		ArrayList<String> columns = getColumns();
		for (String row : rows) {
			for (String col : columns) {
				if (this.isMetaDataColumn(col) || this.isMetaDataRow(row)) {
					String value = getCell(row, col);
					workSheet.addCell(row, col, value);
				} else {
					String value = getCell(row, col);
					try {
						Double d = Double.parseDouble(value);
						if (d == 0.0) {
							d = zeroValue;
						} else {
							d = Math.log(d) / Math.log(base);
						}
						workSheet.addCell(row, col, d + "");
					} catch (Exception e) {
						workSheet.addCell(row, col, value);
					}

				}
			}
		}

		ArrayList<String> metadataRows = this.getMetaDataRows();
		ArrayList<String> metadataColumns = this.getMetaDataColumns();
		workSheet.setMetaDataColumns(metadataColumns);
		workSheet.setMetaDataRows(metadataRows);
		return workSheet;
	}

	/**
	 * Swap the row and columns returning a new worksheet
	 *
	 * @return
	 * @throws Exception
	 */
	public WorkSheet swapRowAndColumns() throws Exception {

		WorkSheet swappedWorkSheet = new WorkSheet(getColumns(), getRows());
		for (String row : getRows()) {
			for (String col : getColumns()) {
				String value = getCell(row, col);
				swappedWorkSheet.addCell(col, row, value);
			}
		}

		ArrayList<String> metadataRows = this.getMetaDataRows();
		ArrayList<String> metadataColumns = this.getMetaDataColumns();
		swappedWorkSheet.setMetaDataColumns(metadataRows);
		swappedWorkSheet.setMetaDataRows(metadataColumns);
		return swappedWorkSheet;
	}

	static CompactCharSequence[][] getAllValuesCompactCharSequence(File fileName, char delimiter) throws Exception {
		FileInputStream fi = new FileInputStream(fileName);
		return getAllValuesCompactCharSequence(fi, delimiter);
	}

	/**
	 * All support for loading from a jar file
	 *
	 * @param is
	 * @param delimiter
	 * @return
	 * @throws Exception
	 */
	static CompactCharSequence[][] getAllValuesCompactCharSequence(InputStream is, char delimiter) throws Exception {
		// FileReader reader = new FileReader(fileName);

		BufferedReader br = new BufferedReader(new InputStreamReader(is));


		ArrayList<CompactCharSequence[]> rows = new ArrayList<CompactCharSequence[]>();

		String line = br.readLine();
		int numcolumns = -1;
		while (line != null) {
			String[] d = line.split(String.valueOf(delimiter));
			if (numcolumns == -1) {
				numcolumns = d.length;
			}
			CompactCharSequence[] ccs = new CompactCharSequence[d.length];
			for (int i = 0; i < d.length; i++) {
				ccs[i] = new CompactCharSequence(d[i]);
			}
			rows.add(ccs);

			line = br.readLine();
		}
		br.close();
		// reader.close();

		CompactCharSequence[][] data = new CompactCharSequence[rows.size()][numcolumns];
		for (int i = 0; i < rows.size(); i++) {
			CompactCharSequence[] row = rows.get(i);
			for (int j = 0; j < row.length; j++) { //
				if (row[j].length() > 1 && row[j].charAt(0) == '"') {
					// System.out.println(row[j]);
					if (row[j].length() > 2) {
						row[j] = new CompactCharSequence(row[j].subSequence(1, row[j].length() - 1).toString());
					} else {
						row[j] = new CompactCharSequence("");
					}
				}
				if (j < row.length && j < data[0].length) {
					data[i][j] = row[j];
				}
			}


		}

		return data;

	}

	static String[][] getAllValues(String fileName, char delimiter) throws Exception {
		FileReader reader = new FileReader(fileName);
		BufferedReader br = new BufferedReader(reader);
		ArrayList<String[]> rows = new ArrayList<String[]>();

		String line = br.readLine();
		int numcolumns = -1;
		while (line != null) {
			String[] d = line.split(String.valueOf(delimiter));
			if (numcolumns == -1) {
				numcolumns = d.length;
			}
			rows.add(d);

			line = br.readLine();
		}
		br.close();
		reader.close();

		String[][] data = new String[rows.size()][numcolumns];
		for (int i = 0; i < rows.size(); i++) {
			String[] row = rows.get(i);
			for (int j = 0; j < row.length; j++) {
				if (row[j].startsWith("\"") && row[j].endsWith("\"")) {
					// System.out.println(row[j]);
					row[j] = row[j].substring(1, row[j].length() - 1);
				}
				data[i][j] = row[j];
			}


		}

		return data;

	}

	/**
	 * Combine two work sheets where you join based on rows. Rows that are found
	 * in one but not the other are removed. If the second sheet is meta data
	 * then a meta data column will be added between the two joined columns
	 *
	 * @param w1FileName
	 * @param w2FileName
	 * @param delimitter
	 * @param secondSheetMetaData
	 * @return
	 * @throws Exception
	 */
	static public WorkSheet unionWorkSheetsRowJoin(String w1FileName, String w2FileName, char delimitter, boolean secondSheetMetaData) throws Exception {
		WorkSheet w1 = WorkSheet.readCSV(w1FileName, delimitter);
		WorkSheet w2 = WorkSheet.readCSV(w2FileName, delimitter);
		return unionWorkSheetsRowJoin(w1, w2, secondSheetMetaData);

	}

	/**
	 * * Combine two work sheets where you join based on rows. Rows that are
	 * found in one but not the other are removed. If the second sheet is meta
	 * data then a meta data column will be added between the two joined columns
	 *
	 * @param w1
	 * @param w2
	 * @param secondSheetMetaData
	 * @return
	 * @throws Exception
	 */
	static public WorkSheet unionWorkSheetsRowJoin(WorkSheet w1, WorkSheet w2, boolean secondSheetMetaData) throws Exception {
		ArrayList<String> w1Columns = w1.getColumns();
		ArrayList<String> w2Columns = w2.getColumns();
		ArrayList<String> w1DataColumns = w1.getDataColumns();
		ArrayList<String> w2DataColumns = w2.getDataColumns();
		ArrayList<String> w1MetaDataColumns = w1.getMetaDataColumns();
		ArrayList<String> w2MetaDataColumns = w2.getMetaDataColumns();


		if (secondSheetMetaData) {
			if (!w1.getColumns().contains("META_DATA")) {
				w1DataColumns.add("META_DATA");
			}
		}

		ArrayList<String> joinedColumns = new ArrayList<String>();
		joinedColumns.addAll(w1DataColumns);
		joinedColumns.addAll(w2DataColumns);
		if (!joinedColumns.contains("META_DATA") && (w1MetaDataColumns.size() > 0 || w2MetaDataColumns.size() > 0)) {
			joinedColumns.add("META_DATA");
		}
		for (String column : w1MetaDataColumns) {
			if (!joinedColumns.contains(column)) {
				joinedColumns.add(column);
			}
		}
		for (String column : w2MetaDataColumns) {
			if (!joinedColumns.contains(column)) {
				joinedColumns.add(column);
			}
		}
		ArrayList<String> w1Rows = w1.getRows();
		ArrayList<String> w2Rows = w2.getRows();
		ArrayList<String> rows = new ArrayList<String>();

		HashSet<String> w1Key = new HashSet<String>(w1Rows);
		for (String key : w2Rows) {
			if (w1Key.contains(key)) {
				rows.add(key);
			}
		}

		WorkSheet worksheet = new WorkSheet(rows, joinedColumns);

		for (String row : rows) {
			for (String column : w1Columns) {
				if (column.equals("META_DATA")) {
					continue;
				}
				String value = w1.getCell(row, column);
				worksheet.addCell(row, column, value);
			}
		}

		for (String row : rows) {
			for (String column : w2Columns) {
				if (column.equals("META_DATA")) {
					continue;
				}
				String value = w2.getCell(row, column);
				worksheet.addCell(row, column, value);
			}
		}
		worksheet.setMetaDataColumnsAfterColumn();
		worksheet.setMetaDataRowsAfterRow();
		return worksheet;
	}

	/**
	 * Read a CSV/Tab delimitted file where you pass in the delimiter
	 *
	 * @param fileName
	 * @param delimiter
	 * @return
	 * @throws Exception
	 */
	static public WorkSheet readCSV(String fileName, char delimiter) throws Exception {


		return readCSV(new File(fileName), delimiter);
	}

	static public WorkSheet readCSV(File f, char delimiter) throws Exception {


		return readCSV(new FileInputStream(f), delimiter);
	}

	/**
	 * Read a CSV/Tab delimited file where you pass in the delimiter
	 *
	 * @param f
	 * @param delimiter
	 * @return
	 * @throws Exception
	 */
	static public WorkSheet readCSV(InputStream is, char delimiter) throws Exception {


		CompactCharSequence[][] data = getAllValuesCompactCharSequence(is, delimiter);

		WorkSheet workSheet = new WorkSheet(data);
		workSheet.setMetaDataColumnsAfterColumn();
		workSheet.setMetaDataRowsAfterRow();
		return workSheet;
	}

	/**
	 * Save the worksheet as a csv file
	 *
	 * @param fileName
	 * @throws Exception
	 */
	public void saveCSV(String fileName) throws Exception {
		File f = new File(fileName);
		File parentFile = f.getParentFile();
		if (!parentFile.isDirectory()) {
			parentFile.mkdirs();
		}
		FileOutputStream file = new FileOutputStream(fileName);
		BufferedOutputStream bs = new BufferedOutputStream(file);
		save(bs, ',', false);
		bs.close();
		file.close();
	}

	/**
	 *
	 * @param fileName
	 * @throws Exception
	 */
	public void saveTXT(String fileName) throws Exception {
		File f = new File(fileName);
		File parentFile = f.getParentFile();
		if (!parentFile.isDirectory()) {
			parentFile.mkdirs();
		}
		FileOutputStream file = new FileOutputStream(fileName);
		BufferedOutputStream bs = new BufferedOutputStream(file);
		save(bs, '\t', false);
		bs.close();
		file.close();
	}
	private String rowHeader = "REF";

	/**
	 *
	 * @param value
	 */
	public void setRowHeader(String value) {
		rowHeader = value;
	}

	/**
	 * Add columns from a second worksheet to be joined by common row. If the
	 * appended worksheet doesn't contain a row in the master worksheet then
	 * default value of "" is used. Rows in the appended worksheet not found in
	 * the master worksheet are not added.
	 *
	 * @param worksheet
	 * @throws Exception
	 */
	public void appendWorkSheetColumns(WorkSheet worksheet) throws Exception {

		ArrayList<String> newColumns = worksheet.getColumns();

		this.addColumns(newColumns, "");
		ArrayList<String> rows = this.getRows();
		for (String row : rows) {
			for (String col : newColumns) {
				if (worksheet.isValidRow(row)) {
					String value = worksheet.getCell(row, col);
					this.addCell(row, col, value);
				}

			}
		}



	}

	/**
	 * Add rows from a second worksheet to be joined by common column. If the
	 * appended worksheet doesn't contain a column in the master worksheet then
	 * default value of "" is used. Columns in the appended worksheet not found
	 * in the master worksheet are not added.
	 *
	 * @param worksheet
	 * @throws Exception
	 */
	public void appendWorkSheetRows(WorkSheet worksheet) throws Exception {

		ArrayList<String> newRows = worksheet.getRows();

		this.addRows(newRows, "");
		for (String col : this.getColumns()) {
			if (!worksheet.isValidColumn(col)) {
				continue;
			}
			for (String row : newRows) {
				if (worksheet.isValidColumn(col)) {
					String value = worksheet.getCell(row, col);
					this.addCell(row, col, value);
				}

			}
		}

	}

	/**
	 *
	 * @param outputStream
	 * @param delimitter
	 * @param quoteit
	 * @throws Exception
	 */
	public void save(OutputStream outputStream, char delimitter, boolean quoteit) throws Exception {
		outputStream.write(rowHeader.getBytes());
		//String quote = "\"";

		for (String col : getColumns()) {
			outputStream.write(delimitter);
			if (quoteit) {
				outputStream.write('"');
			}
			outputStream.write(col.getBytes());
			if (quoteit) {
				outputStream.write('"');
			}
		}
		outputStream.write("\r\n".getBytes());
		for (String row : getRows()) {
			if (quoteit) {
				outputStream.write('"');
			}
			outputStream.write(row.getBytes());
			if (quoteit) {
				outputStream.write('"');
			}
			for (String col : getColumns()) {
				// try{
				String value = getCell(row, col);
				outputStream.write(delimitter);
				if (!this.isMetaDataColumn(col) && !this.isMetaDataRow(row)) {
					if (value == null || value.length() == 0 || value.equalsIgnoreCase("null")) {
						value = "NaN";
					}
				} else {
					if (value == null || value.length() == 0 || value.equalsIgnoreCase("null")) {
						value = "";
					}
				}

				outputStream.write(value.getBytes());
				//  }catch(Exception e){
				//      System.out.println(row + " " + col);
				//  }
			}
			outputStream.write("\r\n".getBytes());
		}
	}

	/**
	 * @return the indexColumnName
	 */
	public String getIndexColumnName() {
		return indexColumnName;
	}

	/**
	 * @param indexColumnName the indexColumnName to set
	 */
	public void setIndexColumnName(String indexColumnName) {
		this.indexColumnName = indexColumnName;
	}

	/**
	 * @return the columnLookup
	 */
	public LinkedHashMap<String, HeaderInfo> getColumnLookup() {
		return columnLookup;
	}

	/**
	 * @return the rowLookup
	 */
	public LinkedHashMap<String, HeaderInfo> getRowLookup() {
		return rowLookup;
	}

	/**
	 * @return the metaDataColumnsHashMap
	 */
	public LinkedHashMap<String, String> getMetaDataColumnsHashMap() {
		return metaDataColumnsHashMap;
	}

	/**
	 * @return the metaDataRowsHashMap
	 */
	public LinkedHashMap<String, String> getMetaDataRowsHashMap() {
		return metaDataRowsHashMap;
	}

	/**
	 * @return the rowHeader
	 */
	public String getRowHeader() {
		return rowHeader;
	}
}
