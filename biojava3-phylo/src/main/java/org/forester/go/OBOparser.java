// $Id: OBOparser.java,v 1.15 2009/11/10 23:09:38 cmzmasek Exp $
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

package org.forester.go;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.forester.util.ForesterUtil;

public class OBOparser {

    private final File       _input_file;    ;
    private final ReturnType _return_type;
    private int              _go_term_count;

    public OBOparser( final File input_file, final ReturnType return_type ) {
        switch ( return_type ) {
            case BASIC_GO_TERM:
                break;
            default:
                throw new IllegalArgumentException( "unknown return type: " + return_type );
        }
        _input_file = input_file;
        _return_type = return_type;
        init();
    }

    private GoTerm createNewBasicGoTerm( final String id,
                                         final String name,
                                         final String namespace,
                                         final String is_obsolete,
                                         final String comment,
                                         final String definition,
                                         final Set<String> alt_ids,
                                         final List<GoXRef> go_xrefs,
                                         final List<GoId> super_go_ids,
                                         final List<GoRelationship> go_relationships,
                                         final List<GoSubset> go_subsets ) {
        final GoTerm gt = new BasicGoTerm( id, name, namespace, is_obsolete.trim().toLowerCase().equals( "true" ) );
        ( ( BasicGoTerm ) gt ).setComment( comment );
        ( ( BasicGoTerm ) gt ).setDefinition( definition );
        for( final GoXRef x : go_xrefs ) {
            gt.getGoXRefs().add( x );
        }
        for( final GoId s : super_go_ids ) {
            gt.getSuperGoIds().add( s );
        }
        for( final GoRelationship r : go_relationships ) {
            gt.getGoRelationships().add( r );
        }
        for( final GoSubset sub : go_subsets ) {
            gt.getGoSubsets().add( sub );
        }
        for( final String alt_id : alt_ids ) {
            gt.getAltIds().add( new GoId( alt_id ) );
        }
        ++_go_term_count;
        return gt;
    }

    private void createNewGoTerm( final List<GoTerm> go_terms,
                                  final String id,
                                  final String name,
                                  final String namespace,
                                  final String is_obsolete,
                                  final String comment,
                                  final String definition,
                                  final Set<String> alt_ids,
                                  final List<GoXRef> go_xrefs,
                                  final List<GoId> super_go_ids,
                                  final List<GoRelationship> go_relationships,
                                  final List<GoSubset> go_subsets ) {
        GoTerm gt;
        switch ( getReturnType() ) {
            case BASIC_GO_TERM:
                gt = createNewBasicGoTerm( id,
                                           name,
                                           namespace,
                                           is_obsolete,
                                           comment,
                                           definition,
                                           alt_ids,
                                           go_xrefs,
                                           super_go_ids,
                                           go_relationships,
                                           go_subsets );
                break;
            default:
                throw new AssertionError( "unknown return type: " + getReturnType() );
        }
        go_terms.add( gt );
    }

    public int getGoTermCount() {
        return _go_term_count;
    }

    private File getInputFile() {
        return _input_file;
    }

    private ReturnType getReturnType() {
        return _return_type;
    }

    private void init() {
        setGoTermCount( 0 );
    }

    public List<GoTerm> parse() throws IOException {
        final String error = ForesterUtil.isReadableFile( getInputFile() );
        if ( !ForesterUtil.isEmpty( error ) ) {
            throw new IOException( error );
        }
        final BufferedReader br = new BufferedReader( new FileReader( getInputFile() ) );
        String line;
        final List<GoTerm> go_terms = new ArrayList<GoTerm>();
        int line_number = 0;
        boolean in_term = false;
        String id = "";
        String name = "";
        String namespace = "";
        String def = "";
        String comment = "";
        String is_obsolete = "";
        HashSet<String> alt_ids = new HashSet<String>();
        List<GoId> super_go_ids = new ArrayList<GoId>();
        List<GoXRef> go_xrefs = new ArrayList<GoXRef>();
        List<GoRelationship> go_relationships = new ArrayList<GoRelationship>();
        List<GoSubset> go_subsets = new ArrayList<GoSubset>();
        try {
            while ( ( line = br.readLine() ) != null ) {
                line_number++;
                line = line.trim();
                if ( line.length() < 1 ) {
                    if ( in_term ) {
                        in_term = false;
                    }
                }
                else if ( line.startsWith( "[Term]" ) ) {
                    in_term = true;
                    if ( id.length() > 0 ) {
                        createNewGoTerm( go_terms,
                                         id,
                                         name,
                                         namespace,
                                         is_obsolete,
                                         comment,
                                         def,
                                         alt_ids,
                                         go_xrefs,
                                         super_go_ids,
                                         go_relationships,
                                         go_subsets );
                    }
                    id = "";
                    name = "";
                    namespace = "";
                    alt_ids = new HashSet<String>();
                    def = "";
                    comment = "";
                    is_obsolete = "";
                    super_go_ids = new ArrayList<GoId>();
                    go_xrefs = new ArrayList<GoXRef>();
                    go_relationships = new ArrayList<GoRelationship>();
                    go_subsets = new ArrayList<GoSubset>();
                }
                else if ( in_term && line.startsWith( "id:" ) ) {
                    id = line.substring( 3 ).trim();
                }
                else if ( in_term && line.startsWith( "name:" ) ) {
                    name = line.substring( 5 ).trim();
                }
                else if ( in_term && line.startsWith( "namespace:" ) ) {
                    namespace = line.substring( 10 ).trim();
                }
                else if ( in_term && line.startsWith( "alt_id:" ) ) {
                    alt_ids.add( line.substring( 7 ).trim() );
                }
                else if ( in_term && line.startsWith( "def:" ) ) {
                    def = line.substring( 4 ).trim();
                }
                else if ( in_term && line.startsWith( "is_obsolete:" ) ) {
                    is_obsolete = line.substring( 12 ).trim();
                }
                else if ( in_term && line.startsWith( "comment:" ) ) {
                    comment = line.substring( 8 ).trim();
                }
                else if ( in_term && line.startsWith( "xref:" ) ) {
                    final String s = trimOffComment( line.substring( 5 ).trim() );
                    go_xrefs.add( new BasicGoXRef( s ) );
                }
                else if ( in_term && line.startsWith( "is_a:" ) ) {
                    final String s = trimOffComment( line.substring( 5 ).trim() );
                    super_go_ids.add( new GoId( s ) );
                }
                else if ( in_term && line.startsWith( "relationship:" ) ) {
                    final String s = trimOffComment( line.substring( 13 ).trim() );
                    go_relationships.add( new BasicGoRelationship( s ) );
                }
                else if ( in_term && line.startsWith( "subset:" ) ) {
                    final String s = line.substring( 8 ).trim();
                    go_subsets.add( new BasicGoSubset( s ) );
                }
            } // while ( ( line = br.readLine() ) != null )
        }
        catch ( final Exception e ) {
            throw new IOException( "parsing problem: " + e.getMessage() + " [at line " + line_number + "]" );
        }
        if ( id.length() > 0 ) {
            createNewGoTerm( go_terms,
                             id,
                             name,
                             namespace,
                             is_obsolete,
                             comment,
                             def,
                             alt_ids,
                             go_xrefs,
                             super_go_ids,
                             go_relationships,
                             go_subsets );
        }
        return go_terms;
    }

    private void setGoTermCount( final int go_term_count ) {
        _go_term_count = go_term_count;
    }

    private String trimOffComment( String xref ) {
        final int i = xref.indexOf( '!' );
        if ( i > 0 ) {
            xref = xref.substring( 0, xref.indexOf( '!' ) ).trim();
        }
        return xref;
    }

    public static enum ReturnType {
        BASIC_GO_TERM
    }
}
