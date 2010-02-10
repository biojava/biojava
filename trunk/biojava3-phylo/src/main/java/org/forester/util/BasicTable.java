// $Id: BasicTable.java,v 1.15 2009/11/11 01:36:50 cmzmasek Exp $
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
// WWW: www.phylosoft.org

package org.forester.util;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class BasicTable<E> {

    private Map<String, Map<String, E>> _rows;
    private int                         _max_row;
    private int                         _max_col;

    public BasicTable() {
        init();
    }

    public Map<String, E> getColumnsAsMap( final int key_col, final int value_col ) throws IllegalArgumentException {
        final Map<String, E> map = new HashMap<String, E>();
        for( int row = 0; row < getNumberOfRows(); ++row ) {
            final String key = ( String ) getValue( key_col, row );
            final E value = getValue( value_col, row );
            if ( ( key != null ) && ( value != null ) ) {
                if ( map.containsKey( key ) ) {
                    throw new IllegalArgumentException( "attempt to use non-unique table value as key [" + key + "]" );
                }
                map.put( key, value );
            }
        }
        return map;
    }

    public Map<String, Double> getColumnsAsMapDouble( final int key_col, final int value_col )
            throws IllegalArgumentException, IOException {
        final Map<String, Double> map = new HashMap<String, Double>();
        for( int row = 0; row < getNumberOfRows(); ++row ) {
            final String key = ( String ) getValue( key_col, row );
            double value = 0;
            try {
                value = Double.parseDouble( getValueAsString( value_col, row ) );
            }
            catch ( final NumberFormatException e ) {
                throw new IOException( e );
            }
            if ( key != null ) {
                if ( map.containsKey( key ) ) {
                    throw new IllegalArgumentException( "attempt to use non-unique table value as key [" + key + "]" );
                }
                map.put( key, value );
            }
        }
        return map;
    }

    public int getNumberOfColumns() {
        return _max_col + 1;
    }

    public int getNumberOfRows() {
        return _max_row + 1;
    }

    private Map<String, E> getRow( final int row ) {
        return getRows().get( "" + row );
    }

    private Map<String, Map<String, E>> getRows() {
        return _rows;
    }

    public E getValue( final int col, final int row ) throws IllegalArgumentException {
        if ( ( row > getNumberOfRows() - 1 ) || ( row < 0 ) ) {
            throw new IllegalArgumentException( "value for row (" + row + ") is out of range [number of rows: "
                    + getNumberOfRows() + "]" );
        }
        else if ( ( col >= getNumberOfColumns() ) || ( row < 0 ) ) {
            throw new IllegalArgumentException( "value for column (" + col + ") is out of range [number of columns: "
                    + getNumberOfColumns() + "]" );
        }
        final Map<String, E> row_map = getRow( row );
        if ( ( row_map == null ) || ( row_map.size() < 1 ) ) {
            return null;
        }
        return row_map.get( "" + col );
    }

    public String getValueAsString( final int col, final int row ) throws IllegalArgumentException {
        if ( getValue( col, row ) != null ) {
            return getValue( col, row ).toString();
        }
        return null;
    }

    private void init() {
        _rows = new HashMap<String, Map<String, E>>();
        setMaxCol( -1 );
        setMaxRow( -1 );
    }

    public boolean isEmpty() {
        return getNumberOfRows() <= 0;
    }

    private void setMaxCol( final int max_col ) {
        _max_col = max_col;
    }

    private void setMaxRow( final int max_row ) {
        _max_row = max_row;
    }

    public void setValue( final int col, final int row, final E value ) {
        if ( ( row < 0 ) || ( col < 0 ) ) {
            throw new IllegalArgumentException( "attempt to use negative values for row or column" );
        }
        if ( row > getNumberOfRows() - 1 ) {
            setMaxRow( row );
        }
        if ( col > getNumberOfColumns() - 1 ) {
            setMaxCol( col );
        }
        final String row_key = "" + row;
        Map<String, E> row_map = null;
        if ( getRows().containsKey( row_key ) ) {
            row_map = getRows().get( row_key );
        }
        else {
            row_map = new HashMap<String, E>();
            getRows().put( row_key, row_map );
        }
        row_map.put( "" + col, value );
    }

    @Override
    public String toString() {
        final StringBuffer sb = new StringBuffer();
        for( int row = 0; row < getNumberOfRows(); ++row ) {
            for( int col = 0; col < getNumberOfColumns(); ++col ) {
                sb.append( getValue( col, row ) );
                if ( col < getNumberOfColumns() - 1 ) {
                    sb.append( " " );
                }
            }
            if ( row < getNumberOfRows() - 1 ) {
                sb.append( ForesterUtil.LINE_SEPARATOR );
            }
        }
        return sb.toString();
    }
}
