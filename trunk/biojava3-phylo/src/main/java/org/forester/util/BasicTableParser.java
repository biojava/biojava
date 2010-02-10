// $Id: BasicTableParser.java,v 1.14 2009/04/21 01:44:00 cmzmasek Exp $
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

package org.forester.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

public class BasicTableParser {

    private final static String START_OF_COMMENT_LINE_DEFAULT = "#";

    private BasicTableParser() {
    }

    public static BasicTable<String> parse( final Object source, final String column_delimiter ) throws IOException {
        return BasicTableParser.parse( source, column_delimiter, false, START_OF_COMMENT_LINE_DEFAULT, false ).get( 0 );
    }

    public static BasicTable<String> parse( final Object source,
                                            final String column_delimiter,
                                            final boolean use_first_separator_only ) throws IOException {
        return BasicTableParser.parse( source,
                                       column_delimiter,
                                       use_first_separator_only,
                                       START_OF_COMMENT_LINE_DEFAULT,
                                       false ).get( 0 );
    }

    public static List<BasicTable<String>> parse( final Object source,
                                                  final String column_delimiter,
                                                  final boolean use_first_separator_only,
                                                  final String start_of_comment_line,
                                                  final boolean tables_separated_by_single_string_line )
            throws IOException {
        final BufferedReader reader = ForesterUtil.obtainReader( source );
        final List<BasicTable<String>> tables = new ArrayList<BasicTable<String>>();
        BasicTable<String> table = new BasicTable<String>();
        int row = 0;
        String line;
        boolean saw_first_table = false;
        final boolean use_start_of_comment_line = !( ForesterUtil.isEmpty( start_of_comment_line ) );
        while ( ( line = reader.readLine() ) != null ) {
            line = line.trim();
            if ( saw_first_table
                    && ( ForesterUtil.isEmpty( line ) || ( tables_separated_by_single_string_line && ( line
                            .indexOf( column_delimiter ) < 0 ) ) ) ) {
                if ( !table.isEmpty() ) {
                    tables.add( table );
                }
                table = new BasicTable<String>();
                row = 0;
            }
            else if ( !ForesterUtil.isEmpty( line )
                    && ( !use_start_of_comment_line || !line.startsWith( start_of_comment_line ) ) ) {
                saw_first_table = true;
                final StringTokenizer st = new StringTokenizer( line, column_delimiter );
                int col = 0;
                if ( st.hasMoreTokens() ) {
                    table.setValue( col++, row, st.nextToken().trim() );
                }
                if ( !use_first_separator_only ) {
                    while ( st.hasMoreTokens() ) {
                        table.setValue( col++, row, st.nextToken().trim() );
                    }
                }
                else {
                    final StringBuffer rest = new StringBuffer();
                    while ( st.hasMoreTokens() ) {
                        rest.append( st.nextToken() );
                    }
                    table.setValue( col++, row, rest.toString().trim() );
                }
                ++row;
            }
        }
        if ( !table.isEmpty() ) {
            tables.add( table );
        }
        reader.close();
        return tables;
    }
}
