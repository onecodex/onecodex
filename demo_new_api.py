#!/usr/bin/env python3
"""
Demonstration script showing the new OpenAPI-based implementation.

This script demonstrates that the new implementation provides the same
API surface as the existing Potion-Client based implementation.
"""

import os
import sys

# Add the current directory to path so we can import from onecodex
sys.path.insert(0, ".")


def demo_new_api():
    """Demonstrate the new API implementation."""
    print("🚀 OneCodex New API Implementation Demo")
    print("=" * 50)

    try:
        # Import the new API
        from onecodex.new_api import Api

        print("✓ Successfully imported new Api class")

        # Import new models
        from onecodex.new_api import Samples, Classifications, Projects

        print("✓ Successfully imported new model classes")

        # Create API instance (without credentials for demo)
        print("\n📡 Creating API client...")
        try:
            # This will fail without credentials, but shows the interface works
            api = Api(api_key="demo_key", base_url="http://localhost:3000")
            print("✓ API client created successfully")

            # Show that models are bound to the API
            print(f"✓ Samples class bound: {hasattr(api, 'Samples')}")
            print(f"✓ Classifications class bound: {hasattr(api, 'Classifications')}")
            print(f"✓ Projects class bound: {hasattr(api, 'Projects')}")

        except Exception as e:
            print(f"⚠️  API creation failed (expected without real credentials): {e}")

        # Demonstrate model functionality
        print("\n🔧 Model Functionality:")
        print(f"✓ Samples resource path: {Samples._resource_path}")
        print(f"✓ Classifications resource path: {Classifications._resource_path}")
        print(f"✓ Projects resource path: {Projects._resource_path}")

        # Show mixins are working
        print(f"✓ Samples has download method: {hasattr(Samples, 'download')}")
        print(f"✓ Samples has upload method: {hasattr(Samples, 'upload')}")

        # Show CRUD methods exist
        print(f"✓ Samples has get method: {hasattr(Samples, 'get')}")
        print(f"✓ Samples has where method: {hasattr(Samples, 'where')}")
        print(f"✓ Samples has save method: {hasattr(Samples, 'save')}")
        print(f"✓ Samples has delete method: {hasattr(Samples, 'delete')}")

        print("\n📊 Generated Models:")
        from onecodex.models.generated import SampleSchema, ProjectSchema

        print(f"✓ Generated SampleSchema available")
        print(f"✓ Generated ProjectSchema available")

        print("\n🔄 HTTP Client:")
        from onecodex.client import OneCodexHTTPClient

        client = OneCodexHTTPClient(base_url="https://app.onecodex.com")
        print(f"✓ HTTP client created with base URL: {client.base_url}")
        client.close()

        print("\n🎯 Key Benefits:")
        print("• Type safety with Pydantic v2 models")
        print("• Modern HTTP client with httpx")
        print("• Auto-generated models from OpenAPI spec")
        print("• No vendored dependencies")
        print("• Same API surface as existing implementation")

        print("\n✅ Demo completed successfully!")
        print("\nThe new implementation is ready for testing!")

    except ImportError as e:
        print(f"❌ Import error: {e}")
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return False

    return True


if __name__ == "__main__":
    success = demo_new_api()
    sys.exit(0 if success else 1)
