// $Id: BasicMsa.java,v 1.9 2009/10/26 23:29:40 cmzmasek Exp $
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

package org.forester.development;

import java.util.ArrayList;
import java.util.Arrays;

public class BasicMsa implements Msa {

    private final ArrayList<Sequence> _sequences;

    public BasicMsa() {
        _sequences = new ArrayList<Sequence>();
    }

    public void addSequence( final Sequence sequence, final boolean check_for_equal_length ) {
        if ( check_for_equal_length && ( getNumberOfSequences() > 0 ) && ( getLength() != sequence.getLength() ) ) {
            final String msg = "Attempt to a sequence of length " + sequence.getLength() + " to alignment of length "
                    + getLength();
            throw new IllegalArgumentException( msg );
        }
        getSequences().add( sequence );
    }

    public Msa copy() {
        return createSubAlignment( null, null );
    }

    // set sequences to null or empty to get all sequences
    // set positions to null or empty to get all posistions
    public Msa createSubAlignment( final int[] sequences, final int[] positions ) {
        final Msa msa = new BasicMsa();
        boolean return_all_seqs = false;
        boolean return_all_pos = false;
        if ( ( sequences == null ) || ( sequences.length < 1 ) ) {
            return_all_seqs = true;
        }
        else {
            Arrays.sort( sequences );
        }
        if ( ( positions == null ) || ( positions.length < 1 ) ) {
            return_all_pos = true;
        }
        else {
            Arrays.sort( positions );
        }
        for( int s = 0; s < getNumberOfSequences(); ++s ) {
            if ( return_all_seqs || ( Arrays.binarySearch( sequences, s ) >= 0 ) ) {
                for( int p = 0; p < getLength(); ++p ) {
                    if ( return_all_pos || ( Arrays.binarySearch( positions, p ) >= 0 ) ) {
                        // FIXME finish me
                    }
                }
            }
        }
        return msa;
    }

    public int[] findSequence( final String name, final boolean case_sensitive, final boolean partial_match ) {
        // TODO Auto-generated method stub
        return null;
    }

    public int getLength() {
        if ( isEmpty() ) {
            return 0;
        }
        return ( getSequence( 0 ) ).getLength();
    }

    public int getNumberOfSequences() {
        return getSequences().size();
    }

    @Override
    public char getResidueAt( final int col, final int row ) {
        // TODO Auto-generated method stub
        return 0;
    }

    public Sequence getSequence( final int i ) {
        return getSequences().get( i );
    }

    private ArrayList<Sequence> getSequences() {
        return _sequences;
    }

    public boolean isAppearsAligned() {
        if ( !isEmpty() ) {
            final int l = getSequence( 0 ).getLength();
            for( int i = 0; i < getNumberOfSequences(); ++i ) {
                if ( l != getSequence( i ).getLength() ) {
                    return false;
                }
            }
        }
        return true;
    }

    public boolean isEmpty() {
        return getNumberOfSequences() < 1;
    }

    public void removeSequence( final int i ) {
        getSequences().remove( i );
    }
}
