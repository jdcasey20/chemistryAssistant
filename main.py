from dataclasses import dataclass, field
from drafter import *
import periodictable


"""The periodictable library provides atomic masses and parses
   chemical formulas. It is used in this program to calculate molar masses and
   count atoms when checking whether equations are balanced.
   This avoids manually storing atomic masses for every element.
"""



# -----------------------------------------
# DATA CLASSES
# -----------------------------------------

def calculate_molar_mass(formula: str) -> float:
    """Return the molar mass in g/mol of a chemical formula using periodictable.
    """
    chemical_formula = periodictable.formula(formula)
    return chemical_formula.mass


def calculate_moles(mass: float, formula: str) -> float:
    """Convert a mass in grams to moles using the chemical formula.
    """
    molar_mass = calculate_molar_mass(formula)
    return mass / molar_mass



def get_coefficient(equation: str, chemical: str) -> int:
    """Find the stoichiometric coefficient of a substance in an equation.
    
    Return zero if the equation or substance cannot be matched.
    """
    equation = equation.replace(" ", "")
    equation = equation.replace("→", "->")

    sides = equation.split("->")

    if len(sides) != 2:
        return 0

    substances = sides[0].split("+") + sides[1].split("+")

    for substance in substances:
        index = 0

        while index < len(substance) and substance[index].isdigit():
            index += 1

        coefficient_text = substance[:index]
        formula = substance[index:]

        if formula == chemical:
            if coefficient_text == "":
                return 1
            else:
                return int(coefficient_text)

    return 0



def find_limiting_reactant(equation: str,
                           reactant_1: str, mass_1: float,
                           reactant_2: str, mass_2: float,
                           reactant_3: str = "",
                           mass_3: float = 0.0) -> str:

    """Identify the limiting reactant from the supplied reactant masses.
    
    Compare available moles divided by stoichiometric coefficients.
    An optional third reactant may be included.
    """
    reactants = [
        (reactant_1, mass_1),
        (reactant_2, mass_2)
    ]

    if reactant_3 != "":
        reactants.append((reactant_3, mass_3))

    # Infinity ensures the first valid ratio becomes the initial minimum.
    smallest_ratio = float("inf")
    limiting = ""

    for formula, mass in reactants:
        coefficient = get_coefficient(equation, formula)

        if coefficient == 0:
            return "Invalid reactant or equation"

        moles = calculate_moles(mass, formula)
        ratio = moles / coefficient

        if ratio < smallest_ratio:
            smallest_ratio = ratio
            limiting = formula

    return limiting

def calculate_theoretical_yield(equation: str,
                                reactant_1: str, mass_1: float,
                                reactant_2: str, mass_2: float,
                                product: str,
                                reactant_3: str = "",
                                mass_3: float = 0.0) -> float:

    """Calculate the maximum product mass in grams from reactant masses.
    
    Use the smallest reactant mole-to-coefficient ratio to determine
    the moles of product that can form.
    """
    reactants = [
        (reactant_1, mass_1),
        (reactant_2, mass_2)
    ]

    if reactant_3 != "":
        reactants.append((reactant_3, mass_3))

    product_coefficient = get_coefficient(equation, product)

    if product_coefficient == 0:
        return 0.0

    # Infinity ensures the first valid ratio becomes the initial minimum.
    smallest_ratio = float("inf")

    for formula, mass in reactants:
        coefficient = get_coefficient(equation, formula)

        if coefficient == 0:
            return 0.0

        moles = calculate_moles(mass, formula)
        ratio = moles / coefficient

        if ratio < smallest_ratio:
            smallest_ratio = ratio

    product_moles = smallest_ratio * product_coefficient
    product_mass = product_moles * calculate_molar_mass(product)

    return product_mass


def calculate_percent_yield(actual_yield: float,
                            theoretical_yield: float) -> float:

    """Return actual yield divided by theoretical yield as a percentage.
    
    Return zero if the theoretical yield is not positive.
    """
    if theoretical_yield <= 0:
        return 0.0

    percent_yield = (actual_yield / theoretical_yield) * 100

    return percent_yield


def calculate_reaction_heat(equation: str,
                            limiting_reactant: str,
                            limiting_mass: float,
                            enthalpy_change: float) -> float:

    """Calculate heat transferred in kJ from the limiting reactant.
    
    The enthalpy change is in kJ per mole of reaction as written.
    """
    coefficient = get_coefficient(equation, limiting_reactant)

    if coefficient == 0:
        return 0.0

    moles = calculate_moles(limiting_mass, limiting_reactant)

    # Enthalpy is given per mole of reaction as the equation is written.
    reaction_amount = moles / coefficient

    heat = reaction_amount * enthalpy_change

    return heat

def validate_inputs(reactant_1_formula: str,
                    reactant_1_amount: str,
                    reactant_2_formula: str,
                    reactant_2_amount: str,
                    reactant_3_formula: str,
                    reactant_3_amount: str,
                    product_name: str,
                    actual_yield: str,
                    calculate_heat: bool,
                    enthalpy_change: str) -> str:

    # Check required chemical formulas
    """Check required inputs and numerical values for a calculation.
    
    Return an explanatory error message, or an empty string if valid.
    """
    if reactant_1_formula.strip() == "":
        return "Please enter the first reactant formula."

    if reactant_2_formula.strip() == "":
        return "Please enter the second reactant formula."

    if product_name.strip() == "":
        return "Please enter a product formula."

    # Check required masses
    if reactant_1_amount.strip() == "":
        return "Please enter the mass of reactant 1."

    if reactant_2_amount.strip() == "":
        return "Please enter the mass of reactant 2."

    # Validate the optional third reactant
    if reactant_3_formula.strip() != "" and reactant_3_amount.strip() == "":
        return "Please enter the mass of reactant 3."

    if reactant_3_formula.strip() == "" and reactant_3_amount.strip() != "":
        return "Please enter a formula for reactant 3."

    # Check that masses are valid numbers
    masses = [reactant_1_amount, reactant_2_amount]

    if reactant_3_formula.strip() != "":
        masses.append(reactant_3_amount)

    for mass in masses:
        try:
            number = float(mass)
        except ValueError:
            return "Reactant masses must be numbers."

        if number <= 0 or number == float("inf") or number != number:
            return "Reactant masses must be positive, finite numbers."

    # Actual yield is optional
    if actual_yield.strip() != "":
        try:
            actual = float(actual_yield)
        except ValueError:
            return "Actual yield must be a number."

        if actual < 0 or actual == float("inf") or actual == float("-inf") or actual != actual:
            return "Actual yield must be a nonnegative, finite number."

    # Enthalpy is required only when heat calculation is selected
    if calculate_heat:
        if enthalpy_change.strip() == "":
            return "Please enter the reaction enthalpy change."

        try:
            delta_h = float(enthalpy_change)
        except ValueError:
            return "Enthalpy change must be a number."

        if delta_h == float("inf") or delta_h == float("-inf") or delta_h != delta_h:
            return "Enthalpy change must be a finite number."

    return ""





def count_atoms(side: str) -> dict[str, float]:
    # Store the total number of atoms of each element
    """Count atoms of each element on one side of a chemical equation.
    
    Account for stoichiometric coefficients and return a dictionary
    mapping element symbols to total atom counts.
    """
    atom_counts = {}

    substances = side.replace(" ", "").split("+")

    for substance in substances:
        index = 0

        # Find the coefficient at the beginning
        while index < len(substance) and substance[index].isdigit():
            index += 1

        if index == 0:
            coefficient = 1
        else:
            coefficient = int(substance[:index])

        if coefficient <= 0:
            raise ValueError("Coefficients must be positive.")

        formula = substance[index:]

        # periodictable converts a chemical formula
        # into a dictionary of elements and atom counts
        chemical = periodictable.formula(formula)

        for element, amount in chemical.atoms.items():
            symbol = str(element)

            if symbol not in atom_counts:
                atom_counts[symbol] = 0

            atom_counts[symbol] += amount * coefficient

    return atom_counts


def validate_equation(equation: str,
                      reactant_1: str,
                      reactant_2: str,
                      reactant_3: str,
                      product: str) -> str:

    """Check equation syntax, atom balance, and supplied chemical formulas.
    
    Verify that entered reactants occur on the left and the product
    occurs on the right. Return an error message or an empty string.
    """
    equation = equation.replace("→", "->")
    equation = equation.replace(" ", "")

    sides = equation.split("->")

    if len(sides) != 2:
        return "Use -> to separate reactants and products."

    left_side = sides[0]
    right_side = sides[1]

    if left_side == "" or right_side == "":
        return "The equation must have reactants and products."

    # Check whether all formulas and coefficients are valid
    try:
        reactant_atoms = count_atoms(left_side)
        product_atoms = count_atoms(right_side)

        # Validate the individual formulas entered by the user
        formulas = [reactant_1, reactant_2, product]

        if reactant_3 != "":
            formulas.append(reactant_3)

        for formula in formulas:
            periodictable.formula(formula)

    except (ValueError, KeyError, TypeError):
        return "One or more chemical formulas are invalid."

    # Verify that the same number of each atom
    # appears on both sides of the equation
    if reactant_atoms != product_atoms:
        return "The chemical equation is not balanced."

    # Check that reactants are on the left
    left_formulas = []

    for substance in left_side.split("+"):
        index = 0

        while index < len(substance) and substance[index].isdigit():
            index += 1

        left_formulas.append(substance[index:])

    if reactant_1 not in left_formulas:
        return "Reactant 1 is not on the left side of the equation."

    if reactant_2 not in left_formulas:
        return "Reactant 2 is not on the left side of the equation."

    if reactant_3 != "" and reactant_3 not in left_formulas:
        return "Reactant 3 is not on the left side of the equation."

    # Check that the product is on the right
    right_formulas = []

    for substance in right_side.split("+"):
        index = 0

        while index < len(substance) and substance[index].isdigit():
            index += 1

        right_formulas.append(substance[index:])

    if product not in right_formulas:
        return "The product is not on the right side of the equation."

    return ""




@dataclass
class ReactionRecord:
    equation: str
    limiting_reactant: str
    theoretical_yield: float
    actual_yield: float
    percent_yield: float
    heat_generated: float
    has_actual_yield: bool = False
    has_heat_calculation: bool = False

"""The State dataclass stores information between website pages.
   default_factory=list creates a separate empty reaction history
   for each new State instance, instead of sharing one list.
"""

@dataclass
class State:
    reaction_history: list[ReactionRecord] = field(default_factory=list)
    
    """ This field can contain either a ReactionRecord or None.
        None means the user has not completed a calculation yet.
    """
    current_reaction: ReactionRecord | None = None
    selected_reaction_index: int = -1

    # Chemical equation information
    equation: str = ""

    # First reactant
    reactant_1_formula: str = ""
    reactant_1_amount: str = ""

    # Second reactant
    reactant_2_formula: str = ""
    reactant_2_amount: str = ""

    # Third reactant
    reactant_3_formula: str = ""
    reactant_3_amount: str = ""

    # Product information
    product_name: str = ""
    actual_yield: str = ""

    # Additional settings
    calculate_heat: bool = False
    enthalpy_change: str = ""
    message: str = ""

# -----------------------------------------
# ROUTE 1: HOME PAGE
# -----------------------------------------

@route
def index(state: State) -> Page:
    """Display the home page and navigation to the calculator and guide.
    """
    return Page(state, [
        Header("CHEMISTRY REACTION CALCULATOR", 3),
        Text("Welcome to the Chemistry Reaction Calculator!"),
        Text("Use this website to analyze chemical reactions "
             "and perform important stoichiometric calculations."),
        Header("WHAT CAN THIS WEBSITE DO?", 2),
        BulletedList([
            "Identify the limiting reactant",
            "Calculate theoretical yield",
            "Calculate percent yield",
            "Calculate heat generated by a reaction",
            "Save and review previous calculations"
        ]),
        Button("Start Calculating", "/reaction_input"),
        Button("Reaction History", "/history"),
        Button("Chemistry Guide", "/guide")
    ])


# -----------------------------------------
# ROUTE 2: REACTION INPUT PAGE
# -----------------------------------------


@route
def reaction_input(state: State) -> Page:
    """Display the reaction entry form with values saved in State.
    """
    return Page(state, [
        Header("ENTER YOUR CHEMICAL REACTION", 1),

        Text("Enter a balanced chemical equation and the quantities of your reactants."),

        Header("Chemical Equation", 2),
        TextBox("equation", state.equation),

        Header("Reactant 1", 2),
        Text("Chemical Formula:"),
        TextBox("reactant_1_formula", state.reactant_1_formula),
        Text("Amount in grams:"),
        TextBox("reactant_1_amount", state.reactant_1_amount, "number"),

        Header("Reactant 2", 2),
        Text("Chemical Formula:"),
        TextBox("reactant_2_formula", state.reactant_2_formula),
        Text("Amount in grams:"),
        TextBox("reactant_2_amount", state.reactant_2_amount, "number"),

        Header("Reactant 3 (Optional)", 2),
        Text("Chemical Formula:"),
        TextBox("reactant_3_formula", state.reactant_3_formula),
        Text("Amount in grams:"),
        TextBox("reactant_3_amount", state.reactant_3_amount, "number"),

        Header("Product Information", 2),
        Text("Product Formula:"),
        TextBox("product_name", state.product_name),
        Text("Actual Yield in grams (optional):"),
        TextBox("actual_yield", state.actual_yield, "number"),

        Text("Calculate Heat Generated:"),
        CheckBox("calculate_heat", state.calculate_heat),
        Text("Reaction Enthalpy Change (kJ per mole of reaction):"),
        TextBox("enthalpy_change", state.enthalpy_change),

        Button("Calculate Reaction", "/calculate"),
        Button("Back Home", "/")
    ])

# -----------------------------------------
# ROUTE 3: CALCULATION PAGE
# -----------------------------------------


@route
def calculate(state: State,
              equation: str,
              reactant_1_formula: str,
              reactant_1_amount: str,
              reactant_2_formula: str,
              reactant_2_amount: str,
              reactant_3_formula: str,
              reactant_3_amount: str,
              product_name: str,
              actual_yield: str,
              calculate_heat: bool,
              enthalpy_change: str) -> Page:

    """Validate reaction inputs, calculate results, and save a record.
    
    Update State with the submitted values and reaction history,
    then display either an error page or a completion page.
    """
    # Save the information entered by the user
    state.equation = equation

    state.reactant_1_formula = reactant_1_formula
    state.reactant_1_amount = reactant_1_amount

    state.reactant_2_formula = reactant_2_formula
    state.reactant_2_amount = reactant_2_amount

    state.reactant_3_formula = reactant_3_formula
    state.reactant_3_amount = reactant_3_amount

    state.product_name = product_name
    state.actual_yield = actual_yield
    state.calculate_heat = calculate_heat
    state.enthalpy_change = enthalpy_change

   
   
    error = validate_inputs(
        reactant_1_formula,
        reactant_1_amount,
        reactant_2_formula,
        reactant_2_amount,
        reactant_3_formula,
        reactant_3_amount,
        product_name,
        actual_yield,
        calculate_heat,
        enthalpy_change
    )

    if error != "":
        state.message = error
        return Page(state, [
            Header("INVALID INPUT", 1),
            Text(error),
            Button("Try Again", "/reaction_input")
        ])
    
    
    equation_error = validate_equation(
        equation,
        reactant_1_formula,
        reactant_2_formula,
        reactant_3_formula,
        product_name
    )

    if equation_error != "":
        state.message = equation_error
        return Page(state, [
            Header("INVALID EQUATION", 1),
            Text(equation_error),
            Button("Try Again", "/reaction_input")
        ])
   
   
       # Convert the entered masses to numbers
    mass_1 = float(reactant_1_amount)
    mass_2 = float(reactant_2_amount)

    mass_3 = 0.0
    if reactant_3_formula != "":
        mass_3 = float(reactant_3_amount)

    # Identify the limiting reactant
    limiting = find_limiting_reactant(
        equation,
        reactant_1_formula, mass_1,
        reactant_2_formula, mass_2,
        reactant_3_formula, mass_3
    )

    # Calculate theoretical yield
    theoretical = calculate_theoretical_yield(
        equation,
        reactant_1_formula, mass_1,
        reactant_2_formula, mass_2,
        product_name,
        reactant_3_formula, mass_3
    )

    # Calculate percent yield if actual yield was entered
    actual = 0.0
    percent = 0.0

    if actual_yield != "":
        actual = float(actual_yield)
        percent = calculate_percent_yield(actual, theoretical)

    # Calculate heat if requested
    heat = 0.0

    if calculate_heat and enthalpy_change != "":
        delta_h = float(enthalpy_change)

        if limiting == reactant_1_formula:
            limiting_mass = mass_1
        elif limiting == reactant_2_formula:
            limiting_mass = mass_2
        else:
            limiting_mass = mass_3

        heat = calculate_reaction_heat(
            equation,
            limiting,
            limiting_mass,
            delta_h
    )

    # Save the completed calculation
    record = ReactionRecord(
        equation,
        limiting,
        theoretical,
        actual,
        percent,
        heat,
        actual_yield != "",
        calculate_heat
    )

    state.current_reaction = record
    state.reaction_history.append(record)

    state.message = "Your reaction was calculated successfully."

    return Page(state, [
        Header("CALCULATION COMPLETE", 1),
        Text(state.message),
        Button("View Results", "/results"),
        Button("Try Another Reaction", "/reaction_input")
    ])


# -----------------------------------------
# ROUTE 4: RESULTS PAGE
# -----------------------------------------


@route
def results(state: State) -> Page:

    """Display the latest calculated reaction and its available results.
    
    Indicate when optional actual-yield or heat values were not calculated.
    """
    if state.current_reaction is None:
        return Page(state, [
            Header("NO RESULTS AVAILABLE", 1),
            Text("Please enter a reaction first."),
            Button("Start Calculating", "/reaction_input")
        ])

    reaction = state.current_reaction
    
    if state.actual_yield == "":
        actual_text = "Not provided"
        percent_text = "Not calculated"
    else:
        actual_text = str(reaction.actual_yield) + " grams"
        percent_text = str(round(reaction.percent_yield, 2)) + "%"

    if state.calculate_heat:
        heat_text = str(round(reaction.heat_generated, 2)) + " kJ"
    else:
        heat_text = "Not calculated"
    
    return Page(state, [
        Header("REACTION RESULTS", 1),

        Header("Chemical Equation", 2),
        Text(state.equation),

        Header("Reactants", 2),
        Text("Reactant 1: " + state.reactant_1_formula),
        Text("Amount: " + state.reactant_1_amount + " grams"),

        Text("Reactant 2: " + state.reactant_2_formula),
        Text("Amount: " + state.reactant_2_amount + " grams"),

        Text("Reactant 3: " + state.reactant_3_formula),
        Text("Amount: " + state.reactant_3_amount + " grams"),

        Header("Product Information", 2),
        Text("Product: " + state.product_name),

        Header("Calculated Results", 2),

        Text("Limiting Reactant: " + reaction.limiting_reactant),

        Text("Theoretical Yield: " +
             str(round(reaction.theoretical_yield, 2)) + " grams"),

        Text("Actual Yield: " + actual_text),
        Text("Percent Yield: " + percent_text),
        Text("Heat Transferred: " + heat_text),

        Button("New Calculation", "/reaction_input"),
        Button("View History", "/history"),
        Button("Home", "/")
    ])



# -----------------------------------------
# ROUTE 5: REACTION HISTORY PAGE
# -----------------------------------------

@route
def history(state: State) -> Page:

    """Display saved reactions and buttons to open their details.
    
    Iterate over the ReactionRecord objects in State.reaction_history.
    """
    history_content = [
        Header("Reaction History"),
        Text("Review your previously calculated reactions.")
    ]

    if len(state.reaction_history) == 0:
        history_content.append(
            Text("You have not completed any calculations yet.")
        )
    else:
        for index, reaction in enumerate(state.reaction_history):
            history_content.append(
                Header("Reaction " + str(index + 1), 2)
            )

            history_content.append(
                Text(reaction.equation)
            )

          
            # Pass this record's list index to the selection route.
            history_content.append(
                Button("View Details", "select_reaction",
                       [Argument("index", index)])
            )

    history_content.append(Button("Home", "/"))
    history_content.append(Button("New Calculation", "/reaction_input"))

    return Page(state, history_content)


@route
def select_reaction(state: State, index: int) -> Page:
    """Save the selected reaction index and open its details page.
    """
    state.selected_reaction_index = index
    return reaction_details(state)

# -----------------------------------------
# ROUTE 6: REACTION DETAILS PAGE
# -----------------------------------------

@route
def reaction_details(state: State) -> Page:

    """Display the selected saved reaction, or a not-found message.
    
    Show optional results only when they were supplied or calculated.
    """
    index = state.selected_reaction_index

    if 0 <= index < len(state.reaction_history):
        reaction = state.reaction_history[index]
        if reaction.has_actual_yield:
            actual_text = str(reaction.actual_yield) + " grams"
            percent_text = str(round(reaction.percent_yield, 2)) + "%"
        else:
            actual_text = "Not provided"
            percent_text = "Not calculated"

        if reaction.has_heat_calculation:
            heat_text = str(round(reaction.heat_generated, 2)) + " kJ"
        else:
            heat_text = "Not calculated"

        return Page(state, [
            Header("Reaction Details", 2),
            Text("Equation: " + reaction.equation),
            Text("Limiting Reactant: " + reaction.limiting_reactant),
            Text("Theoretical Yield: " + str(round(reaction.theoretical_yield, 2)) + " grams"),
            Text("Actual Yield: " + actual_text),
            Text("Percent Yield: " + percent_text),
            Text("Heat Transferred: " + heat_text),
            Button("Back to History", "/history")
        ])

    else:
        return Page(state, [
            Header("Reaction Not Found", 3),
            Text("The requested reaction does not exist."),
            Button("Back to History", "/history")
        ])


# -----------------------------------------
# ROUTE 7: CHEMISTRY GUIDE PAGE
# -----------------------------------------

@route
def guide(state: State) -> Page:
    """Display explanations of the chemistry calculations on the website.
    """
    return Page(state, [
        Header("Chemistry Guide"),

        Header("Limiting Reactant"),
        Text("The limiting reactant is the reactant that "
             "gets completely consumed first in a chemical reaction."),

        Header("Theoretical Yield"),
        Text("The theoretical yield is the maximum amount "
             "of product that can be produced based on the limiting reactant."),

        Header("Percent Yield"),
        Text("Percent yield compares the actual amount of product "
             "obtained to the theoretical amount."),
        Text("Percent Yield = (Actual Yield / Theoretical Yield) x 100"),

        Header("Heat Generated"),
        Text("Heat generated can be calculated using the enthalpy "
             "change of a reaction and the amount of reaction occurring."),

        Button("Back Home", "/")
    ])


# -----------------------------------------
# START THE WEBSITE
# -----------------------------------------

def run_unit_tests() -> None:
    """Test chemistry helpers and each Drafter route using isolated State objects.

    These tests run before the website starts. They verify returned pages,
    state changes, numerical calculations, and validation error paths.
    """
    equation = "2H2 + O2 -> 2H2O"

    # Chemistry helper tests, including boundary and invalid cases.
    assert abs(calculate_molar_mass("H2O") - 18.015) < 0.1
    assert abs(calculate_moles(18.015, "H2O") - 1.0) < 0.01
    assert get_coefficient(equation, "H2") == 2
    assert get_coefficient(equation, "O2") == 1
    assert get_coefficient(equation, "CO2") == 0
    assert find_limiting_reactant(equation, "H2", 10, "O2", 20) == "O2"
    assert find_limiting_reactant(equation, "H2", 1, "O2", 20) == "H2"
    yield_mass = calculate_theoretical_yield(
        equation, "H2", 10, "O2", 20, "H2O"
    )
    assert 22 < yield_mass < 23
    assert calculate_theoretical_yield(
        equation, "H2", 10, "O2", 20, "CO2"
    ) == 0
    assert calculate_percent_yield(8, 20) == 40
    assert calculate_percent_yield(8, 0) == 0
    heat = calculate_reaction_heat(equation, "O2", 20, -571.6)
    assert -360 < heat < -350
    assert validate_inputs(
        "H2", "10", "O2", "20", "", "", "H2O", "", False, ""
    ) == ""
    assert validate_inputs(
        "H2", "-1", "O2", "20", "", "", "H2O", "", False, ""
    ) != ""
    assert validate_inputs(
        "H2", "10", "O2", "20", "", "", "H2O", "", True, ""
    ) != ""
    assert validate_inputs(
        "H2", "10", "O2", "20", "", "", "H2O", "-inf", False, ""
    ) != ""
    assert count_atoms("2H2+O2") == count_atoms("2H2O")
    assert validate_equation(equation, "H2", "O2", "", "H2O") == ""
    assert validate_equation(
        "H2 + O2 -> H2O", "H2", "O2", "", "H2O"
    ) != ""

    # Test all eight routes with fresh State instances to avoid test overlap.
    assert index(State()) is not None
    assert reaction_input(State()) is not None
    assert results(State()) is not None  # No-current-reaction branch.
    assert history(State()) is not None  # Empty-history branch.
    assert reaction_details(State()) is not None  # Invalid-index branch.
    assert guide(State()) is not None

    state = State()
    # Route 3: valid submission creates a saved record.
    assert calculate(
        state, equation, "H2", "10", "O2", "20", "", "", "H2O",
        "8", True, "-571.6"
    ) is not None
    assert len(state.reaction_history) == 1
    assert state.current_reaction is not None
    assert state.current_reaction.limiting_reactant == "O2"
    assert state.current_reaction.has_actual_yield
    assert state.current_reaction.has_heat_calculation
    assert results(state) is not None  # Populated-results branch.
    assert history(state) is not None  # Populated-history branch.
    assert select_reaction(state, 0) is not None  # Route 6: select.
    assert state.selected_reaction_index == 0
    assert reaction_details(state) is not None  # Valid-record branch.

    # Route 3: invalid submission must not add a record.
    assert calculate(
        state, equation, "H2", "-1", "O2", "20", "", "", "H2O",
        "", False, ""
    ) is not None
    assert len(state.reaction_history) == 1
    assert state.message == "Reactant masses must be positive, finite numbers."

    # Route 3: an unbalanced equation must not add a record.
    assert calculate(
        state, "H2 + O2 -> H2O", "H2", "10", "O2", "20", "", "",
        "H2O", "", False, ""
    ) is not None
    assert len(state.reaction_history) == 1

    # Optional fields must be stored distinctly from numerical zero.
    optional_state = State()
    assert calculate(
        optional_state, equation, "H2", "10", "O2", "20", "", "",
        "H2O", "", False, ""
    ) is not None
    assert not optional_state.current_reaction.has_actual_yield
    assert not optional_state.current_reaction.has_heat_calculation
    assert results(optional_state) is not None
    assert select_reaction(optional_state, 0) is not None


# Tests execute before the server starts; any failure gives a line number.
run_unit_tests()

set_website_title("Chemistry Assistant")
set_site_information(
    author="James Casey",
    description= "Analyze chemical reactions, calculate yields and heat, identify limiting reactants, and save results.",
    sources="Thanks to the Drafter docs and my study group.",
    planning=["plan.pdf", "design.pdf"],
    links=["https://github.com/jdcasey20/chemistryAssistant"]
)
hide_debug_information()
set_website_framed(False)
start_server(State())
