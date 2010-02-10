// $Id: NexusCharactersParser.java,v 1.6 2009/10/26 23:29:40 cmzmasek Exp $
//
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
// WWW: www.phylosoft.org/

package org.forester.io.parsers.nexus;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.forester.io.parsers.ParserUtils;
import org.forester.io.parsers.PhylogenyParserException;
import org.forester.util.ForesterUtil;

public class NexusCharactersParser {

    final private static String charstatelabels = NexusConstants.CHARSTATELABELS.toLowerCase();
    private Object              _nexus_source;
    private String[]            _char_state_labels;

    public String[] getCharStateLabels() {
        return _char_state_labels;
    }

    private Object getNexusSource() {
        return _nexus_source;
    }

    public void parse() throws IOException {
        reset();
        final BufferedReader reader = ParserUtils.createReader( getNexusSource() );
        String line;
        boolean in_charstatelabels = false;
        final List<String> labels_list = new ArrayList<String>();
        int counter = 1;
        while ( ( line = reader.readLine() ) != null ) {
            line = line.trim();
            if ( ( line.length() > 0 ) && !line.startsWith( "#" ) && !line.startsWith( ">" ) ) {
                if ( line.toLowerCase().startsWith( charstatelabels ) ) {
                    in_charstatelabels = true;
                }
                else if ( in_charstatelabels ) {
                    String label = line;
                    if ( label.indexOf( ' ' ) > 0 ) {
                        final String[] s = label.split( "\\s+" );
                        label = s[ 1 ];
                        int count = -1;
                        try {
                            count = Integer.parseInt( s[ 0 ] );
                        }
                        catch ( final NumberFormatException ex ) {
                            throw new NexusFormatException( "failed to parse character label number from: " + line );
                        }
                        if ( count != counter ) {
                            throw new NexusFormatException( "character label numbers are not in order, current line: "
                                    + line );
                        }
                    }
                    ++counter;
                    label = label.replaceAll( "[\\s;\"',]+", "" );
                    if ( !ForesterUtil.isEmpty( label ) ) {
                        if ( labels_list.contains( label ) ) {
                            throw new NexusFormatException( "character label [" + label + "] is not unique" );
                        }
                        labels_list.add( label );
                    }
                }
                if ( line.endsWith( ";" ) ) {
                    in_charstatelabels = false;
                }
            }
        }
        setCharStateLabels( new String[ labels_list.size() ] );
        int i = 0;
        for( final String label : labels_list ) {
            getCharStateLabels()[ i++ ] = label;
        }
    }

    private void reset() {
        setCharStateLabels( new String[ 0 ] );
    }

    private void setCharStateLabels( final String[] char_state_labels ) {
        _char_state_labels = char_state_labels;
    }

    public void setSource( final Object nexus_source ) throws PhylogenyParserException, IOException {
        if ( nexus_source == null ) {
            throw new PhylogenyParserException( getClass() + ": attempt to parse null object." );
        }
        _nexus_source = nexus_source;
    }
}
