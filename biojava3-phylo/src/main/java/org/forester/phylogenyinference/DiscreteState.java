// $Id: DiscreteState.java,v 1.2 2008/03/09 00:11:14 cmzmasek Exp $
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
// WWW: www.phylosoft.org/forester

package org.forester.phylogenyinference;

import org.forester.util.ForesterUtil;

public class DiscreteState implements Comparable<DiscreteState> {

    final private String _id;

    public DiscreteState( final String id ) {
        if ( ForesterUtil.isEmpty( id ) ) {
            throw new IllegalArgumentException( "attempt to create new discrete state from empty or null string" );
        }
        _id = id.trim();
    }

    @Override
    public int compareTo( final DiscreteState discrete_state ) {
        if ( this == discrete_state ) {
            return 0;
        }
        return getId().compareTo( discrete_state.getId() );
    }

    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( o == null ) {
            throw new IllegalArgumentException( "attempt to check discrete state equality to null" );
        }
        else if ( o.getClass() != this.getClass() ) {
            throw new IllegalArgumentException( "attempt to check discrete stateequality to " + o + " [" + o.getClass()
                    + "]" );
        }
        else {
            return getId().equals( ( ( DiscreteState ) o ).getId() );
        }
    }

    public String getId() {
        return _id;
    }

    @Override
    public int hashCode() {
        return getId().hashCode();
    }

    @Override
    public String toString() {
        return getId();
    }
}