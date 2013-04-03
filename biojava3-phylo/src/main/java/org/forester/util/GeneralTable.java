// $Id: GeneralTable.java,v 1.5 2009/10/26 23:29:40 cmzmasek Exp $
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

import java.text.NumberFormat;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

public class GeneralTable<IDENTIFIER_TYPE, VALUE_TYPE> {

    private Map<IDENTIFIER_TYPE, Map<IDENTIFIER_TYPE, VALUE_TYPE>> _rows;
    private SortedSet<IDENTIFIER_TYPE>                             _row_identifiers;
    private SortedSet<IDENTIFIER_TYPE>                             _column_identifiers;

    public GeneralTable() {
        init();
    }

    public SortedSet<IDENTIFIER_TYPE> getColumnIdentifiers() {
        return _column_identifiers;
    }

    private Map<IDENTIFIER_TYPE, VALUE_TYPE> getRow( final IDENTIFIER_TYPE row ) {
        return getRows().get( row );
    }

    public SortedSet<IDENTIFIER_TYPE> getRowIdentifiers() {
        return _row_identifiers;
    }

    private Map<IDENTIFIER_TYPE, Map<IDENTIFIER_TYPE, VALUE_TYPE>> getRows() {
        return _rows;
    }

    public VALUE_TYPE getValue( final IDENTIFIER_TYPE col, final IDENTIFIER_TYPE row ) throws IllegalArgumentException {
        final Map<IDENTIFIER_TYPE, VALUE_TYPE> row_map = getRow( row );
        if ( ( row_map == null ) || ( row_map.size() < 1 ) ) {
            return null;
        }
        return row_map.get( col );
    }

    public String getValueAsString( final IDENTIFIER_TYPE col, final IDENTIFIER_TYPE row )
            throws IllegalArgumentException {
        final VALUE_TYPE value = getValue( col, row );
        return ( value == null ? "" : getValue( col, row ).toString() );
    }

    private void init() {
        _rows = new HashMap<IDENTIFIER_TYPE, Map<IDENTIFIER_TYPE, VALUE_TYPE>>();
        _row_identifiers = new TreeSet<IDENTIFIER_TYPE>();
        _column_identifiers = new TreeSet<IDENTIFIER_TYPE>();
    }

    public void setValue( final IDENTIFIER_TYPE col, final IDENTIFIER_TYPE row, final VALUE_TYPE value ) {
        getColumnIdentifiers().add( col );
        getRowIdentifiers().add( row );
        Map<IDENTIFIER_TYPE, VALUE_TYPE> row_map = null;
        if ( getRows().containsKey( row ) ) {
            row_map = getRows().get( row );
        }
        else {
            row_map = new HashMap<IDENTIFIER_TYPE, VALUE_TYPE>();
            getRows().put( row, row_map );
        }
        row_map.put( col, value );
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        sb.append( "\t" );
        for( final IDENTIFIER_TYPE col : getColumnIdentifiers() ) {
            sb.append( col.toString() );
            sb.append( "\t" );
        }
        sb.append( ForesterUtil.LINE_SEPARATOR );
        for( final IDENTIFIER_TYPE row : getRowIdentifiers() ) {
            sb.append( row.toString() );
            sb.append( "\t" );
            for( final IDENTIFIER_TYPE col : getColumnIdentifiers() ) {
                sb.append( getValueAsString( col, row ) );
                sb.append( "\t" );
            }
            sb.append( ForesterUtil.LINE_SEPARATOR );
        }
        return sb.toString();
    }

    public String toString( final NumberFormat number_format ) {
        final StringBuilder sb = new StringBuilder();
        sb.append( "\t" );
        for( final IDENTIFIER_TYPE col : getColumnIdentifiers() ) {
            sb.append( col.toString() );
            sb.append( "\t" );
        }
        sb.append( ForesterUtil.LINE_SEPARATOR );
        for( final IDENTIFIER_TYPE row : getRowIdentifiers() ) {
            sb.append( row.toString() );
            sb.append( "\t" );
            for( final IDENTIFIER_TYPE col : getColumnIdentifiers() ) {
                try {
                    sb.append( number_format.format( getValue( col, row ) ) );
                }
                catch ( final IllegalArgumentException e ) {
                    sb.append( getValueAsString( col, row ) );
                }
                sb.append( "\t" );
            }
            sb.append( ForesterUtil.LINE_SEPARATOR );
        }
        return sb.toString();
    }
}