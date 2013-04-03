// $Id: OntologizerResult.java,v 1.8 2009/01/13 19:49:32 cmzmasek Exp $
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

package org.forester.go.etc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.forester.go.GoId;
import org.forester.util.ForesterUtil;

/*
 * 
 * Note: this class has a natural ordering that is inconsistent with equals.
 */
public class OntologizerResult implements Comparable<OntologizerResult> {

    final private GoId    _goid;
    final private int     _pop_total;
    final private int     _pop_term;
    final private int     _study_total;
    final private int     _study_term;
    final private int     _pop_family;
    final private int     _study_family;
    final private int     _nparents;
    final private boolean _is_trivial;
    final private double  _p;
    final private double  _p_adjusted;
    final private double  _p_min;
    final private TYPE    _type;

    private OntologizerResult( final String s ) {
        if ( ForesterUtil.isEmpty( s ) ) {
            throw new IllegalArgumentException( "result string is null or empty" );
        }
        final String[] tokens = s.split( "\t" );
        if ( ( tokens.length != 9 ) && ( tokens.length != 11 ) && ( tokens.length != 12 ) ) {
            throw new IllegalArgumentException( "result string [" + s + "] has unexpected format" );
        }
        _goid = new GoId( tokens[ 0 ] );
        _pop_total = Integer.parseInt( tokens[ 1 ] );
        _pop_term = Integer.parseInt( tokens[ 2 ] );
        _study_total = Integer.parseInt( tokens[ 3 ] );
        _study_term = Integer.parseInt( tokens[ 4 ] );
        if ( tokens.length == 11 ) {
            // Topology Elim
            // ID Pop.total Pop.term Study.total Study.term Pop.family Study.family is.trivial p p.adjusted p.min
            _type = TYPE.TOPOLOGY;
            _pop_family = Integer.parseInt( tokens[ 5 ] );
            _study_family = Integer.parseInt( tokens[ 6 ] );
            _is_trivial = Boolean.parseBoolean( tokens[ 7 ] );
            _p = Double.parseDouble( tokens[ 8 ] );
            _p_adjusted = Double.parseDouble( tokens[ 9 ] );
            _p_min = Double.parseDouble( tokens[ 10 ] );
            _nparents = -1;
        }
        else if ( tokens.length == 9 ) {
            // Term for Term
            // ID Pop.total Pop.term Study.total Study.term p p.adjusted p.min name
            _type = TYPE.TERM_FOR_TERM;
            _pop_family = -1;
            _study_family = -1;
            _nparents = -1;
            _is_trivial = false;
            _p = Double.parseDouble( tokens[ 5 ] );
            _p_adjusted = Double.parseDouble( tokens[ 6 ] );
            _p_min = Double.parseDouble( tokens[ 7 ] );
        }
        else {
            // Parent Child Union
            // ID Pop.total Pop.term Study.total Study.term Pop.family Study.family nparents is.trivial p p.adjusted p.min
            _type = TYPE.PARENT_CHILD;
            _pop_family = Integer.parseInt( tokens[ 5 ] );
            _study_family = Integer.parseInt( tokens[ 6 ] );
            _nparents = Integer.parseInt( tokens[ 7 ] );
            _is_trivial = Boolean.parseBoolean( tokens[ 8 ] );
            _p = Double.parseDouble( tokens[ 9 ] );
            _p_adjusted = Double.parseDouble( tokens[ 10 ] );
            _p_min = Double.parseDouble( tokens[ 11 ] );
        }
    }

    @Override
    public int compareTo( final OntologizerResult o ) {
        if ( this == o ) {
            return 0;
        }
        else if ( getPAdjusted() < o.getPAdjusted() ) {
            return -1;
        }
        else if ( getPAdjusted() > o.getPAdjusted() ) {
            return 1;
        }
        else {
            return 0;
        }
    }

    public GoId getGoId() {
        return _goid;
    }

    public int getNParents() {
        return _nparents;
    }

    public double getP() {
        return _p;
    }

    public double getPAdjusted() {
        return _p_adjusted;
    }

    public double getPMin() {
        return _p_min;
    }

    public int getPopFamily() {
        return _pop_family;
    }

    public int getPopTerm() {
        return _pop_term;
    }

    public int getPopTotal() {
        return _pop_total;
    }

    public int getStudyFamily() {
        return _study_family;
    }

    public int getStudyTerm() {
        return _study_term;
    }

    public int getStudyTotal() {
        return _study_total;
    }

    public TYPE getType() {
        return _type;
    }

    public boolean isTrivial() {
        return _is_trivial;
    }

    public static List<OntologizerResult> parse( final File input ) throws IOException {
        final String error = ForesterUtil.isReadableFile( input );
        if ( !ForesterUtil.isEmpty( error ) ) {
            throw new IOException( error );
        }
        final BufferedReader br = new BufferedReader( new FileReader( input ) );
        String line;
        final List<OntologizerResult> results = new ArrayList<OntologizerResult>();
        int line_number = 0;
        try {
            while ( ( line = br.readLine() ) != null ) {
                line_number++;
                line = line.trim();
                if ( line.startsWith( "GO:" ) ) {
                    results.add( new OntologizerResult( line ) );
                }
            }
        }
        catch ( final Exception e ) {
            throw new IOException( "parsing problem [at line " + line_number + "] in [" + input + "]: "
                    + e.getMessage() );
        }
        return results;
    }

    public static enum TYPE {
        TOPOLOGY, TERM_FOR_TERM, PARENT_CHILD;
    }
}
