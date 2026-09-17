include: "workflow.smk"

rule report:
    """Render the final pipeline report"""
    input:
        rules.report__workflow__render.output,
