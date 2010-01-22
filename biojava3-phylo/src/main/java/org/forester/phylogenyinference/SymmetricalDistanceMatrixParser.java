// $Id: SymmetricalDistanceMatrixParser.java,v 1.4 2007/12/04 09:18:23 cmzmasek
// Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.phylogenyinference;

import java.io.IOException;
import java.util.List;

import org.forester.util.BasicTable;
import org.forester.util.BasicTableParser;
import org.forester.util.ForesterUtil;

/*
 * This can read full, lower triangular, and upper triangular distance matrices.
 * In the case of a full matrix, the lower triangular values are used. Format
 * (by example): id1 0 id2 0.3 0 id3 0.4 0.4 0
 * 
 * OR
 * 
 * id1 id2 0.3 id3 0.4 0.4
 * 
 * Numbers before are after the data are ignored.
 * 
 * 
 * 
 * 
 * @author Christian M Zmasek
 */
public class SymmetricalDistanceMatrixParser implements DistanceMatrixParser {

    private final static InputMatrixType INPUT_MATRIX_TYPE_DEFAULT = InputMatrixType.LOWER_TRIANGLE;
    private final static String          COMMENT                   = "#";
    private final static String          VALUE_SEPARATOR           = " ";
    private int                          _matrix_size;
    private InputMatrixType              _input_matrix_type;

    private SymmetricalDistanceMatrixParser() {
        init();
    }

    private void checkValueIsZero( final BasicTable<String> table, final int row, final int i, final int start_row )
            throws IOException {
        double d = 0.0;
        final String table_value = table.getValue( i, row + start_row );
        if ( ForesterUtil.isEmpty( table_value ) ) {
            throw new IOException( "value is null or empty at [" + ( i - 1 ) + ", " + row + "]" );
        }
        try {
            d = Double.parseDouble( table_value );
        }
        catch ( final NumberFormatException e ) {
            throw new IOException( "illegal format for distance [" + table_value + "] at [" + ( i - 1 ) + ", " + row
                    + "]" );
        }
        if ( !ForesterUtil.isEqual( 0.0, d ) ) {
            throw new IOException( "attempt to use non-zero diagonal value [" + table_value + "] at [" + ( i - 1 )
                    + ", " + row + "]" );
        }
    }

    private InputMatrixType getInputMatrixType() {
        return _input_matrix_type;
    }

    private int getMatrixSize() {
        return _matrix_size;
    }

    private void init() {
        setInputMatrixType( INPUT_MATRIX_TYPE_DEFAULT );
        reset();
    }

    public DistanceMatrix[] parse( final Object source ) throws IOException {
        reset();
        final List<BasicTable<String>> tables = BasicTableParser.parse( source, VALUE_SEPARATOR, false, COMMENT, true );
        final DistanceMatrix[] distance_matrices = new DistanceMatrix[ tables.size() ];
        int i = 0;
        for( final BasicTable<String> table : tables ) {
            distance_matrices[ i++ ] = transform( table );
        }
        return distance_matrices;
    }

    private void reset() {
        setMatrixSize( -1 );
    }

    public void setInputMatrixType( final InputMatrixType input_matrix_type ) {
        _input_matrix_type = input_matrix_type;
    }

    private void setMatrixSize( final int matrix_size ) {
        _matrix_size = matrix_size;
    }

    private void transferValue( final BasicTable<String> table,
                                final DistanceMatrix distance_matrix,
                                final int row,
                                final int col,
                                final int start_row,
                                final int col_offset ) throws IOException {
        double d = 0.0;
        final String table_value = table.getValue( col, row + start_row );
        if ( ForesterUtil.isEmpty( table_value ) ) {
            throw new IOException( "value is null or empty at [" + ( col - 1 ) + ", " + row + "]" );
        }
        try {
            d = Double.parseDouble( table_value );
        }
        catch ( final NumberFormatException e ) {
            throw new IOException( "illegal format for distance [" + table_value + "] at [" + ( col - 1 ) + ", " + row
                    + "]" );
        }
        distance_matrix.setValue( col - 1 + col_offset, row, d );
    }

    private DistanceMatrix transform( final BasicTable<String> table ) throws IllegalArgumentException, IOException {
        boolean first_line_is_size = false;
        if ( table.getNumberOfColumns() < 3 ) {
            throw new IllegalArgumentException( "attempt to create distance matrix with with less than 3 columns [columns: "
                    + table.getNumberOfColumns() + ", rows: " + table.getNumberOfRows() + "]" );
        }
        if ( table.getNumberOfColumns() == table.getNumberOfRows() ) {
            first_line_is_size = true;
        }
        else if ( table.getNumberOfColumns() != table.getNumberOfRows() + 1 ) {
            throw new IllegalArgumentException( "attempt to create distance matrix with illegal dimensions [columns: "
                    + table.getNumberOfColumns() + ", rows: " + table.getNumberOfRows() + "]" );
        }
        final DistanceMatrix distance_matrix = new BasicSymmetricalDistanceMatrix( table.getNumberOfColumns() - 1 );
        int start_row = 0;
        if ( first_line_is_size ) {
            start_row = 1;
        }
        for( int row = 0; row < table.getNumberOfRows() - start_row; row++ ) {
            distance_matrix.setIdentifier( row, table.getValue( 0, row + start_row ) );
            switch ( getInputMatrixType() ) {
                case LOWER_TRIANGLE:
                    for( int col = 1; col <= row; ++col ) {
                        transferValue( table, distance_matrix, row, col, start_row, 0 );
                    }
                    checkValueIsZero( table, row, row + 1, start_row );
                    break;
                case UPPER_TRIANGLE:
                    for( int col = 1; col < ( table.getNumberOfColumns() - row ); ++col ) {
                        transferValue( table, distance_matrix, row, col, start_row, row );
                    }
                    break;
                default:
                    throw new AssertionError( "unkwnown input matrix type [" + getInputMatrixType() + "]" );
            }
        }
        if ( getMatrixSize() < 1 ) {
            setMatrixSize( distance_matrix.getSize() );
        }
        else if ( getMatrixSize() != distance_matrix.getSize() ) {
            throw new IOException( "attempt to use matrices of unequal size: [" + getMatrixSize() + "] vs ["
                    + distance_matrix.getSize() + "]" );
        }
        return distance_matrix;
    }

    public static SymmetricalDistanceMatrixParser createInstance() {
        return new SymmetricalDistanceMatrixParser();
    }

    public enum InputMatrixType {
        UPPER_TRIANGLE, LOWER_TRIANGLE
    }
}
