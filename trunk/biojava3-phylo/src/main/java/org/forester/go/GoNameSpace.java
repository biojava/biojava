// $Id: GoNameSpace.java,v 1.15 2009/04/28 01:59:44 cmzmasek Exp $
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

public class GoNameSpace {

    public final String           MOLECULAR_FUNCTION_STR = "molecular_function";
    public final String           BIOLOGICAL_PROCESS_STR = "biological_process";
    public final String           CELLULAR_COMPONENT_STR = "cellular_component";
    public final String           UNASSIGNED_STR         = "unassigned";
    private final GoNamespaceType _type;

    public GoNameSpace( final GoNamespaceType type ) {
        _type = type;
    };

    public GoNameSpace( final String type ) {
        if ( type.toLowerCase().equals( MOLECULAR_FUNCTION_STR ) ) {
            _type = GoNamespaceType.MOLECULAR_FUNCTION;
        }
        else if ( type.toLowerCase().equals( BIOLOGICAL_PROCESS_STR ) ) {
            _type = GoNamespaceType.BIOLOGICAL_PROCESS;
        }
        else if ( type.toLowerCase().equals( CELLULAR_COMPONENT_STR ) ) {
            _type = GoNamespaceType.CELLULAR_COMPONENT;
        }
        else if ( type.toLowerCase().equals( UNASSIGNED_STR ) ) {
            _type = GoNamespaceType.UNASSIGNED;
        }
        else {
            throw new IllegalArgumentException( "unknown GO namespace: " + type );
        }
    }

    @Override
    public boolean equals( final Object o ) {
        if ( this == o ) {
            return true;
        }
        else if ( ( o == null ) || ( o.getClass() != this.getClass() ) ) {
            return false;
        }
        else {
            return getType() == ( ( GoNameSpace ) o ).getType();
        }
    }

    public GoNamespaceType getType() {
        return _type;
    }

    public boolean isBiologicalProcess() {
        return getType() == GoNamespaceType.BIOLOGICAL_PROCESS;
    }

    public boolean isCellularComponent() {
        return getType() == GoNamespaceType.CELLULAR_COMPONENT;
    }

    public boolean isMolecularFunction() {
        return getType() == GoNamespaceType.MOLECULAR_FUNCTION;
    }

    public boolean isUnassigned() {
        return getType() == GoNamespaceType.UNASSIGNED;
    }

    public String toShortString() {
        switch ( _type ) {
            case BIOLOGICAL_PROCESS:
                return ( "B" );
            case CELLULAR_COMPONENT:
                return ( "C" );
            case MOLECULAR_FUNCTION:
                return ( "M" );
            case UNASSIGNED:
                return ( "?" );
            default:
                throw new IllegalStateException();
        }
    }

    @Override
    public String toString() {
        switch ( _type ) {
            case BIOLOGICAL_PROCESS:
                return ( BIOLOGICAL_PROCESS_STR );
            case CELLULAR_COMPONENT:
                return ( CELLULAR_COMPONENT_STR );
            case MOLECULAR_FUNCTION:
                return ( MOLECULAR_FUNCTION_STR );
            case UNASSIGNED:
                return ( UNASSIGNED_STR );
            default:
                throw new IllegalStateException();
        }
    }

    public static GoNameSpace createBiologicalProcess() {
        return new GoNameSpace( GoNamespaceType.BIOLOGICAL_PROCESS );
    }

    public static GoNameSpace createCellularComponent() {
        return new GoNameSpace( GoNamespaceType.CELLULAR_COMPONENT );
    }

    public static GoNameSpace createMolecularFunction() {
        return new GoNameSpace( GoNamespaceType.MOLECULAR_FUNCTION );
    }

    public static GoNameSpace createUnassigned() {
        return new GoNameSpace( GoNamespaceType.UNASSIGNED );
    }

    public static enum GoNamespaceType {
        MOLECULAR_FUNCTION, BIOLOGICAL_PROCESS, CELLULAR_COMPONENT, UNASSIGNED;
    }
}
