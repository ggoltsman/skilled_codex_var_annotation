from __future__ import annotations

import anyio
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client


async def main() -> None:
    params = StdioServerParameters(
        command=".venv/bin/python",
        args=["-m", "variant_annotation.server"],
    )
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            tools = await session.list_tools()
            print([tool.name for tool in tools.tools])
            query_result = await session.call_tool(
                "query_variant",
                {"chrom": "1", "pos": 11008, "ref": "C", "alt": "G"},
            )
            print(query_result.content[0].text)
            gnomad_result = await session.call_tool(
                "query_gnomad",
                {"chrom": "1", "pos": 11008, "ref": "C", "alt": "G"},
            )
            print(gnomad_result.content[0].text)
            omim_result = await session.call_tool("query_omim", {"gene": "OR4F5"})
            print(omim_result.content[0].text)
            normalize_result = await session.call_tool(
                "normalize_vcf",
                {"input_vcf": "examples/sample.vcf"},
            )
            print(normalize_result.content[0].text)
            summary_result = await session.call_tool(
                "summarize_annotations",
                {"input_vcf": "examples/sample.vcf"},
            )
            print(summary_result.content[0].text)
            annotation_result = await session.call_tool(
                "annotate_vcf",
                {"input_vcf": "examples/sample.vcf"},
            )
            print(annotation_result.content[0].text)


if __name__ == "__main__":
    anyio.run(main)
