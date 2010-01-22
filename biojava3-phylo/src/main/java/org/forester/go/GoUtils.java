// $Id: GoUtils.java,v 1.21 2009/12/02 03:09:08 cmzmasek Exp $
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
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.forester.util.ForesterUtil;

public final class GoUtils {

    private GoUtils() {
    }

    /**
     * This is for counting the how many times each GO term in 'categories'
     * is a (direct or indirect) super term of the GO terms in 'experiment_set'. 
     * 
     * 
     * @param categories the set of super terms to be counted
     * @param experiment_set the list of GO terms to be analyzed
     * @param all_go_terms all terms in the ontology
     * @return
     */
    public static LinkedHashMap<GoId, Integer> countCategories( final List<GoTerm> categories,
                                                                final List<GoTerm> experiment_set,
                                                                final Map<GoId, GoTerm> all_go_terms ) {
        final LinkedHashMap<GoId, Integer> counts = new LinkedHashMap<GoId, Integer>();
        for( final GoTerm experiment_term : experiment_set ) {
            final Set<GoTerm> super_terms = getAllSuperGoTerms( experiment_term.getGoId(), all_go_terms );
            super_terms.add( experiment_term );
            for( final GoTerm cat : categories ) {
                if ( !counts.containsKey( cat.getGoId() ) ) {
                    counts.put( cat.getGoId(), 0 );
                }
                if ( super_terms.contains( cat ) ) {
                    counts.put( cat.getGoId(), 1 + counts.get( cat.getGoId() ) );
                }
            }
        }
        return counts;
    }

    public static LinkedHashMap<GoId, Integer> countCategoriesId( final List<GoId> categories,
                                                                  final List<GoId> experiment_set,
                                                                  final Map<GoId, GoTerm> all_go_terms ) {
        final LinkedHashMap<GoId, Integer> counts = new LinkedHashMap<GoId, Integer>();
        for( final GoId experiment_id : experiment_set ) {
            final Set<GoId> super_ids = new HashSet<GoId>();
            for( final GoTerm term : getAllSuperGoTerms( experiment_id, all_go_terms ) ) {
                super_ids.add( term.getGoId() );
            }
            super_ids.add( experiment_id );
            for( final GoId cat : categories ) {
                if ( !counts.containsKey( cat ) ) {
                    counts.put( cat, 0 );
                }
                if ( super_ids.contains( cat ) ) {
                    counts.put( cat, 1 + counts.get( cat ) );
                }
            }
        }
        return counts;
    }

    public static Map<GoId, GoTerm> createGoIdToGoTermMap( final List<GoTerm> go_terms ) {
        final Map<GoId, GoTerm> go_id_to_term_map = new HashMap<GoId, GoTerm>();
        for( final GoTerm go_term : go_terms ) {
            go_id_to_term_map.put( go_term.getGoId(), go_term );
            for( final GoId alt_id : go_term.getAltIds() ) {
                go_id_to_term_map.put( alt_id, go_term );
            }
        }
        return go_id_to_term_map;
    }

    public static SortedSet<GoTerm> getAllSuperGoTerms( final GoId go_id, final List<GoTerm> go_terms ) {
        final Map<GoId, GoTerm> goid_to_term_map = GoUtils.createGoIdToGoTermMap( go_terms );
        return getAllSuperGoTerms( go_id, goid_to_term_map );
    }

    public static SortedSet<GoTerm> getAllSuperGoTerms( final GoId go_id, final Map<GoId, GoTerm> goid_to_term_map ) {
        if ( !goid_to_term_map.containsKey( go_id ) ) {
            throw new IllegalArgumentException( "GO id [" + go_id + "] not found in GO id to term map" );
        }
        final GoTerm go_term = goid_to_term_map.get( go_id );
        return getAllSuperGoTerms( go_term, goid_to_term_map );
    }

    public static SortedSet<GoId> getAllSuperGoIds( final GoId go_id, final Map<GoId, GoTerm> goid_to_term_map ) {
        final SortedSet<GoId> ids = new TreeSet<GoId>();
        final SortedSet<GoTerm> terms = GoUtils.getAllSuperGoTerms( go_id, goid_to_term_map );
        for( final GoTerm term : terms ) {
            ids.add( term.getGoId() );
        }
        return ids;
    }

    public static SortedSet<GoTerm> getAllSuperGoTerms( final GoTerm go_term, final Map<GoId, GoTerm> goid_to_term_map ) {
        final SortedSet<GoTerm> supers = new TreeSet<GoTerm>();
        getAllSuperGoTerms( go_term, goid_to_term_map, supers );
        return supers;
    }

    private static void getAllSuperGoTerms( final GoTerm go_term,
                                            final Map<GoId, GoTerm> goid_to_term_map,
                                            final Set<GoTerm> supers ) {
        if ( ( go_term.getSuperGoIds() != null ) && ( go_term.getSuperGoIds().size() > 0 ) ) {
            for( final GoId super_go_id : go_term.getSuperGoIds() ) {
                if ( !goid_to_term_map.containsKey( super_go_id ) ) {
                    throw new IllegalArgumentException( "GO id [" + super_go_id + "] not found in GO id to term map" );
                }
                final GoTerm super_go_term = goid_to_term_map.get( super_go_id );
                supers.add( super_go_term );
                getAllSuperGoTerms( super_go_term, goid_to_term_map, supers );
            }
        }
    }

    public static GoTerm getPenultimateGoTerm( final GoTerm go_term, final Map<GoId, GoTerm> map ) {
        GoTerm my_go_term = go_term;
        GoTerm penultimate = my_go_term;
        while ( ( my_go_term.getSuperGoIds() != null ) && ( my_go_term.getSuperGoIds().size() > 0 ) ) {
            penultimate = my_go_term;
            if ( !map.containsKey( my_go_term.getSuperGoIds().get( 0 ) ) ) {
                throw new IllegalArgumentException( "GO-id [" + my_go_term.getSuperGoIds().get( 0 )
                        + "] not found in map" );
            }
            my_go_term = map.get( my_go_term.getSuperGoIds().get( 0 ) );
        }
        return penultimate;
    }

    public static GoTerm getUltimateGoTerm( final GoTerm go_term, final Map<GoId, GoTerm> map ) {
        GoTerm my_go_term = go_term;
        while ( ( my_go_term.getSuperGoIds() != null ) && ( my_go_term.getSuperGoIds().size() > 0 ) ) {
            if ( !map.containsKey( my_go_term.getSuperGoIds().get( 0 ) ) ) {
                throw new IllegalArgumentException( "GO-id [" + my_go_term.getSuperGoIds().get( 0 )
                        + "] not found in map" );
            }
            my_go_term = map.get( my_go_term.getSuperGoIds().get( 0 ) );
        }
        return my_go_term;
    }

    public static SortedMap<String, List<GoId>> parseGoIds( final Object source,
                                                            final String start_of_comment_line,
                                                            final String start_of_label_line ) throws IOException {
        final Pattern label_matcher = Pattern.compile( start_of_label_line + "\\s*(.+?)" );
        final BufferedReader reader = ForesterUtil.obtainReader( source );
        final SortedMap<String, List<GoId>> results = new TreeMap<String, List<GoId>>();
        String line = "";
        String label = "";
        final boolean use_label = !ForesterUtil.isEmpty( start_of_label_line );
        final boolean use_comment = !ForesterUtil.isEmpty( start_of_comment_line );
        List<GoId> current_list = new ArrayList<GoId>();
        while ( ( line = reader.readLine() ) != null ) {
            line = line.trim();
            if ( ForesterUtil.isEmpty( line ) || ( use_comment && line.startsWith( start_of_comment_line ) ) ) {
                continue;
            }
            else if ( use_label && line.startsWith( start_of_label_line ) ) {
                final Matcher matcher = label_matcher.matcher( line );
                if ( matcher.matches() ) {
                    if ( !ForesterUtil.isEmpty( label ) ) {
                        results.put( label, current_list );
                        current_list = new ArrayList<GoId>();
                    }
                    label = matcher.group( 1 );
                }
            }
            else {
                final String[] s = line.split( "\\s+" );
                final GoId id = new GoId( s[ 0 ] );
                current_list.add( id );
            }
        }
        if ( ForesterUtil.isEmpty( label ) ) {
            label = "";
        }
        results.put( label, current_list );
        reader.close();
        return results;
    }
}
